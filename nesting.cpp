#include "nesting.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <thread>
#include <cstdlib>
#include <stdexcept>
#include <functional>
#include <iostream>

// 为MSVC定义PI常量，因为M_PI可能不可用
#ifndef PI
static constexpr double PI = 3.14159265358979323846;
#endif

namespace nesting {

    using Polygon_2 = geo::Polygon_2;

    // 前向声明，提前声明后面会使用的函数
    static bool all_vertices_inside_circle(const Polygon_with_holes_2& pwh, const geo::FT& R);
    static void get_greedy_initial_solution(Layout& layout);
    static void attempt_place_unplaced_part(Layout& L);
    bool minimize_overlap(Layout& layout, volatile bool* requestQuit);

    // 计算利用率：圆形板材计算完全在圆内的零件面积，矩形板材计算已放置零件的面积
    static double compute_util(const Layout& layout, double /*length*/) {
        if (layout.sheets.empty()) return 0.0;  // 如果没有板材，利用率为0
        if (layout.sheets[0].is_circle()) {  // 如果是圆形板材
            geo::FT R = layout.sheets[0].get_radius();  // 获取半径
            if (R <= geo::FT(0)) return 0.0;  // 半径无效则返回0
            double parts_area = 0.0;  // 初始化零件总面积
            if (!layout.sheet_parts.empty()) {  // 如果有零件数据
                for (size_t i = 0; i < layout.sheet_parts[0].size(); ++i) {  // 遍历所有零件
                    const auto& p = layout.sheet_parts[0][i];  // 获取第i个零件
                    if (all_vertices_inside_circle(p.transformed, R)) {  // 如果零件所有顶点都在圆内
                        parts_area += CGAL::to_double(geo::pwh_area(p.transformed));  // 累加零件面积
                    }
                }
            }
            double Rd = CGAL::to_double(R);  // 将半径转换为double
            double sheet_area = PI * Rd * Rd;  // 计算圆形板材面积
            if (sheet_area <= 0.0) return 0.0;  // 面积无效则返回0
            double util = parts_area / sheet_area;  // 计算利用率
            if (util > 1.0) util = 1.0;  // 限制利用率不超过1.0
            return util;  // 返回利用率
        }
        // 矩形板材的回退计算
        double parts_area = 0.0;  // 初始化零件总面积
        const double FAR_THRESHOLD = -1e5;  // 定义远距离阈值，判断零件是否被放置
        if (!layout.sheet_parts.empty()) {  // 如果有零件数据
            for (size_t i = 0; i < layout.sheet_parts[0].size(); ++i) {  // 遍历所有零件
                const auto& p = layout.sheet_parts[0][i];  // 获取第i个零件
                double tx = p.get_translate_double_x();  // 获取零件x坐标
                double ty = p.get_translate_double_y();  // 获取零件y坐标
                if (tx > FAR_THRESHOLD && ty > FAR_THRESHOLD) {  // 如果零件坐标不在远处（表示已放置）
                    parts_area += CGAL::to_double(geo::pwh_area(p.transformed));  // 累加零件面积
                }
            }
        }
        double height = CGAL::to_double(layout.sheets[0].get_height());  // 获取板材高度
        double length = CGAL::to_double(layout.cur_length);  // 获取当前长度
        if (length <= 0.0 || height <= 0.0) return 0.0;  // 尺寸无效则返回0
        const double EPS = 1e-12;  // 小量避免除零
        double rect_area = length * height + EPS;  // 计算矩形板材面积
        double util = parts_area / rect_area;  // 计算利用率
        if (util > 1.0) util = 1.0;  // 限制利用率不超过1.0
        return util;  // 返回利用率
    }

    // 检查多边形所有顶点是否都在圆内
    static bool all_vertices_inside_circle(const Polygon_with_holes_2& pwh, const geo::FT& R) {
        if (R <= geo::FT(0)) return false;  // 半径无效则返回false
        double Rd = CGAL::to_double(R);  // 将半径转换为double
        double R2 = Rd * Rd;  // 计算半径的平方
        // 检查外边界所有顶点
        for (auto vit = pwh.outer_boundary().vertices_begin(); vit != pwh.outer_boundary().vertices_end(); ++vit) {
            double dx = CGAL::to_double(vit->x()) - Rd;  // 计算顶点x坐标到圆心的距离
            double dy = CGAL::to_double(vit->y()) - Rd;  // 计算顶点y坐标到圆心的距离
            if (dx * dx + dy * dy > R2 + 1e-12) return false;  // 如果距离平方大于半径平方（加上容差），返回false
        }
        // 检查所有洞的顶点
        for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
            for (auto vit = hit->vertices_begin(); vit != hit->vertices_end(); ++vit) {
                double dx = CGAL::to_double(vit->x()) - Rd;  // 计算顶点x坐标到圆心的距离
                double dy = CGAL::to_double(vit->y()) - Rd;  // 计算顶点y坐标到圆心的距离
                if (dx * dx + dy * dy > R2 + 1e-12) return false;  // 如果距离平方大于半径平方，返回false
            }
        }
        return true;  // 所有顶点都在圆内，返回true
    }

    // 计算两个多边形之间的NFP（No-Fit Polygon），使用缓存避免重复计算
    NFPCacheValue& comp_nfp(const Polygon_with_holes_2* poly_A,
        const uint32_t rotation_A,
        const uint32_t allowed_rotation_A,
        const Polygon_with_holes_2* poly_B,
        const uint32_t rotation_B,
        const uint32_t allowed_rotation_B,
        Layout& layout) {
        NFPCacheKey key(poly_A, poly_B, rotation_A, rotation_B);  // 创建缓存键
        auto kv = layout.nfp_cache.find(key);  // 在缓存中查找
        if (kv == layout.nfp_cache.end()) {  // 如果缓存未命中
            Transformation scale(CGAL::SCALING, -1);  // 创建缩放变换矩阵，缩放因子为-1（用于镜像）
            auto rotate_A = geo::get_rotate(rotation_A, allowed_rotation_A);  // 获取多边形A的旋转变换
            auto rotate_B = geo::get_rotate(rotation_B, allowed_rotation_B);  // 获取多边形B的旋转变换

            Polygon_with_holes_2 minus_B;  // 存储B的镜像
            Polygon_with_holes_2 rotated_A;  // 存储A的旋转结果

            // 对B应用旋转和镜像变换
            if (rotate_B) minus_B = geo::transform_polygon_with_holes(scale * (*rotate_B), *poly_B);
            else minus_B = geo::transform_polygon_with_holes(scale, *poly_B);

            // 对A应用旋转变换
            if (rotate_A) rotated_A = geo::transform_polygon_with_holes(*rotate_A, *poly_A);
            else rotated_A = *poly_A;

            // 计算Minkowski和得到NFP
            Polygon_with_holes_2 nfp = CGAL::minkowski_sum_2(rotated_A, minus_B);
            geo::strict_simplify(nfp);  // 简化NFP

            auto bbox = nfp.bbox();  // 计算NFP的包围盒
            NFPCacheValue v;  // 创建缓存值对象
            v.nfp = nfp;  // 存储NFP
            // 存储包围盒坐标
            v.xmin = CGAL::to_double(bbox.xmin());
            v.xmax = CGAL::to_double(bbox.xmax());
            v.ymin = CGAL::to_double(bbox.ymin());
            v.ymax = CGAL::to_double(bbox.ymax());

            auto _kv = layout.nfp_cache.emplace(key, v);  // 将结果插入缓存
            return _kv.first->second;  // 返回缓存值
        }
        return kv->second;  // 缓存命中，直接返回
    }

    // Backwards-compatible lightweight overload: forward to full comp_nfp with default (no) rotations
    NFPCacheValue& comp_nfp(const Polygon_with_holes_2* poly_A,
        const Polygon_with_holes_2* poly_B,
        Layout& layout) {
        // default: no rotation and allowed_rotations set to 0
        return comp_nfp(poly_A, 0u, 0u, poly_B, 0u, 0u, layout);
    }

    // 计算穿透距离（penetration depth），使用缓存
    FT comp_pd(NFPCacheValue& v, double px, double py, Layout& layout) {
        auto& original_nfp = v.nfp;  // 获取NFP
        // 检查点是否在NFP包围盒外
        if (px <= v.xmin || px >= v.xmax || py <= v.ymin || py >= v.ymax) return geo::FT_ZERO;

        hash::PDCacheKey key(&original_nfp, px, py);  // 创建PD缓存键
        FT pd;  // 穿透距离
        layout.pd_count++;  // 增加PD计算计数

        auto iter = layout.pd_cache.find(key);  // 在缓存中查找
        if (iter != layout.pd_cache.end()) pd = iter->second;  // 缓存命中
        else {
            layout.pd_miss++;  // 增加缓存未命中计数
            Point_2 relative_point(px, py);  // 创建点对象
            pd = geo::comp_pd(original_nfp, relative_point);  // 计算穿透距离
            layout.pd_cache.emplace(key, pd);  // 将结果插入缓存
        }
        return pd;  // 返回穿透距离
    }

    // 收缩布局：将超出圆形的零件重新放置到圆形内部
    void shrink(Layout& layout) {
        if (layout.sheets.empty()) return;  // 没有板材则返回
        if (!layout.sheets[0].is_circle()) return;  // 不是圆形板材则返回

        std::vector<size_t> moved_indices;  // 存储被移动的零件索引
        geo::FT R = layout.sheets[0].get_radius();  // 获取板材半径
        double sheetR_d = CGAL::to_double(R);  // 转換為double

        // 遍历所有零件
        for (size_t i = 0; i < layout.poly_num; i++) {
            auto& p = layout.sheet_parts[0][i];  // 获取第i个零件

            // 如果零件有顶点在圆形外
            if (!all_vertices_inside_circle(p.transformed, R)) {
                geo::Polygon_2 ifr_poly;  // 内接可行区域多边形

                // 计算零件的最大半径（从第一个顶点到其他顶点的最大距离）
                auto first_v = p.transformed.outer_boundary().vertices_begin();
                double max_sq_d = 0.0;
                for (auto vtx = p.transformed.outer_boundary().vertices_begin(); vtx != p.transformed.outer_boundary().vertices_end(); ++vtx) {
                    double dx = CGAL::to_double(vtx->x() - first_v->x());
                    double dy = CGAL::to_double(vtx->y() - first_v->y());
                    double d2 = dx * dx + dy * dy;
                    if (d2 > max_sq_d) max_sq_d = d2;
                }

                geo::FT rB = geo::FT(std::sqrt(max_sq_d));  // 零件的最大半径
                geo::FT innerR = R - rB;  // 计算内接圆半径

                // 如果存在内接圆
                if (innerR > geo::FT(0)) {
                    const int SEGMENTS = 128;  // 多边形分段数
                    double cx = sheetR_d;  // 圆心x坐标
                    double cy = sheetR_d;  // 圆心y坐标
                    double rr = CGAL::to_double(innerR);  // 内接圆半径

                    // 生成内接多边形（近似圆形）
                    for (int s = 0; s < SEGMENTS; ++s) {
                        double theta = 2.0 * PI * s / SEGMENTS;
                        double xpt = cx + rr * std::cos(theta);
                        double ypt = cy + rr * std::sin(theta);
                        ifr_poly.push_back(Point_2(xpt, ypt));
                    }

                    // 将多边形平移到零件的第一个顶点位置
                    geo::Polygon_2 translated;
                    for (auto vit = ifr_poly.vertices_begin(); vit != ifr_poly.vertices_end(); ++vit) {
                        double nx = CGAL::to_double(vit->x()) - CGAL::to_double(first_v->x());
                        double ny = CGAL::to_double(vit->y()) - CGAL::to_double(first_v->y());
                        translated.push_back(Point_2(nx, ny));
                    }
                    ifr_poly = translated;
                }

                // 如果内接多边形存在，则在内接多边形内随机选择一个位置
                if (ifr_poly.size() > 0) {
                    auto bbox = ifr_poly.bbox();  // 获取包围盒
                    double ux = (double)std::rand() / (double)RAND_MAX;  // 随机x比例
                    double uy = (double)std::rand() / (double)RAND_MAX;  // 随机y比例
                    double random_x = bbox.x_span() * ux + bbox.xmin();  // 计算随机x坐标
                    double random_y = bbox.y_span() * uy + bbox.ymin();  // 计算随机y坐标
                    p.set_translate(random_x, random_y);  // 设置零件新位置
                }
                else {
                    // 如果没有内接多边形，则将零件放在圆心位置
                    p.set_translate(CGAL::to_double(R), CGAL::to_double(R));
                }
                moved_indices.push_back(i);  // 记录被移动的零件索引
            }
        }

        // 更新被移动零件与其他零件的穿透距离
        for (auto& i : moved_indices) {
            auto& p = layout.sheet_parts[0][i];  // 获取被移动的零件
            auto double_px = CGAL::to_double(p.get_translate_ft_x());  // 零件x坐标
            auto double_py = CGAL::to_double(p.get_translate_ft_y());  // 零件y坐标

            // 计算该零件与其他所有零件的穿透距离
            for (size_t j = 0; j < layout.poly_num; j++) {
                if (i == j) continue;  // 跳过自身

                auto& q = layout.sheet_parts[0][j];  // 获取另一个零件
                auto double_qx = q.get_translate_double_x();  // 另一零件x坐标
                auto double_qy = q.get_translate_double_y();  // 另一零件y坐标

                // 计算NFP
                auto& nfp = comp_nfp(q.base, q.get_rotation(), q.allowed_rotations,
                    p.base, p.get_rotation(), p.allowed_rotations, layout);
                // 计算穿透距离
                auto pd = comp_pd(nfp, double_px - double_qx, double_py - double_qy, layout);
                layout.set_pd(i, j, pd);  // 设置穿透距离
            }
        }

        layout.best_utilization = compute_util(layout, 0.0);  // 更新最佳利用率
    }

    // 获取初始解
    void get_init_solu(Layout& layout) {
        if (layout.sheets.empty() || !layout.sheets[0].is_circle())
            throw std::runtime_error("get_init_solu requires a circular sheet");  // 检查是否为圆形板材

        get_greedy_initial_solution(layout);  // 调用贪婪算法获取初始解

        // 如果all_parts为空，则初始化它
        if (layout.all_parts.empty()) {
            std::vector<Polygon_with_holes_2> parts;
            size_t total_parts = layout.sheet_parts.size() ? layout.sheet_parts[0].size() : 0;
            parts.reserve(total_parts);
            for (size_t i = 0; i < total_parts; ++i)
                parts.push_back(layout.sheet_parts[0][i].transformed);  // 添加所有零件
            layout.init_all_parts(parts);  // 初始化all_parts
        }

        const double FAR = -1e6;  // 远距离阈值
        layout.placed_indices.clear();  // 清空已放置零件索引
        layout.total_parts_area = 0.0;  // 重置零件总面积

        // 遍历所有零件，标记已放置和未放置
        for (size_t i = 0; i < layout.poly_num && i < layout.all_parts.size(); ++i) {
            auto& sh = layout.sheet_parts[0][i];  // 获取零件
            double tx = sh.get_translate_double_x();  // 零件x坐标
            double ty = sh.get_translate_double_y();  // 零件y坐标

            // 如果坐标在远处，标记为未放置
            if (tx <= FAR / 2 || ty <= FAR / 2) layout.mark_part_unplaced(i);
            else {
                geo::FT R = layout.sheets[0].get_radius();  // 获取板材半径
                // 如果零件所有顶点都在圆内，标记为已放置
                if (all_vertices_inside_circle(sh.transformed, R)) layout.mark_part_placed(i);
                else layout.mark_part_unplaced(i);  // 否则标记为未放置
            }
        }

        layout.update_cur_length();  // 更新当前长度
        layout.best_length = layout.cur_length;  // 设置最佳长度
        layout.best_result = layout.sheet_parts;  // 设置最佳结果
        layout.best_utilization = compute_util(layout, 0.0);  // 计算最佳利用率
    }

    // 尝试放置未放置的零件
    static void attempt_place_unplaced_part(Layout& L) {
        if (L.sheets.empty() || !L.sheets[0].is_circle()) return;  // 检查是否为圆形板材
        if (L.all_parts.empty()) return;  // 没有零件则返回

        auto unplaced = L.get_unplaced_parts();  // 获取所有未放置的零件
        if (unplaced.empty()) return;  // 没有未放置零件则返回

        size_t pick = (size_t)(std::rand() % unplaced.size());  // 随机选择一个未放置零件
        auto* mp = unplaced[pick];  // 获取指向该零件的指针
        size_t idx = (size_t)(mp - &L.all_parts[0]);  // 计算零件索引
        if (idx >= L.sheet_parts[0].size()) return;  // 索引越界检查

        auto& sh = L.sheet_parts[0][idx];  // 获取零件引用
        uint32_t rotation_i = sh.get_rotation();  // 获取零件旋转角度
        auto rotate = geo::get_rotate(rotation_i, sh.allowed_rotations);  // 获取旋转变换
        // 应用旋转变换（如果有）
        Polygon_with_holes_2 rotated = rotate ? geo::transform_polygon_with_holes(*rotate, *sh.base) : *sh.base;

        // 计算零件的最大半径
        auto first_v = rotated.outer_boundary().vertices_begin();
        double max_sq_d = 0.0;
        for (auto v = rotated.outer_boundary().vertices_begin(); v != rotated.outer_boundary().vertices_end(); ++v) {
            double dx = CGAL::to_double(v->x() - first_v->x());
            double dy = CGAL::to_double(v->y() - first_v->y());
            double d2 = dx * dx + dy * dy;
            if (d2 > max_sq_d) max_sq_d = d2;
        }

        double rB = std::sqrt(max_sq_d);  // 零件的最大半径
        geo::FT Rft = L.sheets[0].get_radius();  // 获取板材半径
        double R = CGAL::to_double(Rft);
        double innerR = R - rB;  // 计算内接圆半径
        if (innerR <= 0.0) return;  // 如果内接圆半径非正，无法放置

        const int RINGS = 6;  // 圆环数
        const int ANGULAR = 32;  // 角度采样数
        double best_tx = 0.0, best_ty = 0.0;  // 最佳位置
        double best_pd = std::numeric_limits<double>::infinity();  // 最小穿透距离
        bool placed = false;  // 是否成功放置

        // 在多个圆环上采样候选位置
        for (int ring = 0; ring <= RINGS && !placed; ++ring) {
            double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);  // 当前圆环半径
            int samples = (ring == 0) ? 1 : ANGULAR;  // 角度采样数（中心点只采样1次）

            for (int s = 0; s < samples; ++s) {
                double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);  // 角度
                double centre_x = R + ring_r * std::cos(theta);  // 候选位置中心x坐标
                double centre_y = R + ring_r * std::sin(theta);  // 候选位置中心y坐标

                // 计算平移量，使零件的第一个顶点位于候选位置
                double tx = centre_x - CGAL::to_double(first_v->x());
                double ty = centre_y - CGAL::to_double(first_v->y());
                sh.set_translate(tx, ty);  // 设置零件位置

                // 检查零件是否完全在圆内
                if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;

                // 计算与所有已放置零件的总穿透距离
                double total_pd = 0.0;
                for (size_t j = 0; j < L.sheet_parts[0].size(); ++j) {
                    if (j == idx) continue;  // 跳过自身

                    double ox = L.sheet_parts[0][j].get_translate_double_x();  // 其他零件x坐标
                    double oy = L.sheet_parts[0][j].get_translate_double_y();  // 其他零件y坐标
                    // 如果其他零件未放置，跳过
                    if (ox < -1e5 || oy < -1e5) continue;

                    // 计算NFP
                    auto& nfp = comp_nfp(L.sheet_parts[0][j].base, L.sheet_parts[0][j].get_rotation(),
                        L.sheet_parts[0][j].allowed_rotations, sh.base, rotation_i,
                        sh.allowed_rotations, L);
                    // 计算穿透距离
                    FT pd_ft = comp_pd(nfp, tx - ox, ty - oy, L);
                    double pd_d = CGAL::to_double(pd_ft);
                    total_pd += pd_d;  // 累加穿透距离

                    // 如果总穿透距离已超过当前最佳值，提前退出
                    if (total_pd > best_pd) break;
                }

                // 如果无穿透（小于等于容差），则记录为最佳位置并标记为已放置
                if (total_pd <= geo::BIAS) {
                    best_tx = tx;
                    best_ty = ty;
                    best_pd = total_pd;
                    placed = true;
                    break;
                }

                // 如果总穿透距离小于当前最佳值，更新最佳位置
                if (total_pd < best_pd) {
                    best_pd = total_pd;
                    best_tx = tx;
                    best_ty = ty;
                }
            }
        }

        // 如果成功放置
        if (placed) {
            sh.set_translate(best_tx, best_ty);  // 设置零件位置
            L.mark_part_placed(idx);  // 标记为已放置

            double placed_x = CGAL::to_double(sh.get_translate_ft_x());  // 零件x坐标
            double placed_y = CGAL::to_double(sh.get_translate_ft_y());  // 零件y坐标

            // 更新该零件与其他所有零件的穿透距离
            for (size_t j = 0; j < L.sheet_parts[0].size(); ++j) {
                if (j == idx) continue;  // 跳过自身

                double ox = L.sheet_parts[0][j].get_translate_double_x();  // 其他零件x坐标
                double oy = L.sheet_parts[0][j].get_translate_double_y();  // 其他零件y坐标

                // 计算NFP
                auto& nfp = comp_nfp(L.sheet_parts[0][j].base, L.sheet_parts[0][j].get_rotation(),
                    L.sheet_parts[0][j].allowed_rotations, sh.base, sh.get_rotation(),
                    sh.allowed_rotations, L);
                // 计算穿透距离
                FT pd = comp_pd(nfp, placed_x - ox, placed_y - oy, L);
                L.set_pd(idx, j, pd);  // 设置穿透距离
            }
        }
        else {
            // 放置失败，将零件移到远处
            sh.set_translate(-1e6, -1e6);
            L.mark_part_unplaced(idx);  // 标记为未放置
        }
    }

    // 最小化重叠（最小化总穿透距离）
    bool minimize_overlap(Layout& layout, volatile bool* requestQuit) {
        if (layout.sheet_parts.empty() || layout.poly_num == 0) return false;

        // 初始化加权因子或其他必要的预处理
        layout.initialize_miu();

        size_t iterations = 0;
        // 构建索引数组以便打乱顺序
        std::vector<size_t> indices(layout.poly_num);
        for (size_t i = 0; i < layout.poly_num; ++i) indices[i] = i;

        // RNG
        std::mt19937_64 rng((unsigned long)clock() ^ (unsigned long)std::hash<std::thread::id>()(std::this_thread::get_id()));
        auto uniform01 = [&](void) {
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            return dist(rng);
        };

        const int TRIES = 32;

        while (iterations < layout.maxIterations) {
            if (*requestQuit) return false;

            bool any_improved = false;

            // 随机遍历次序
            std::shuffle(indices.begin(), indices.end(), rng);

            for (size_t idx_pos = 0; idx_pos < indices.size(); ++idx_pos) {
                if (*requestQuit) return false;
                size_t i = indices[idx_pos];

                double cur_pd = CGAL::to_double(layout.get_one_polygon_pd(i));
                if (cur_pd <= geo::BIAS) continue;

                auto& sh = layout.sheet_parts[0][i];

                // 当前最佳记录
                double best_pd = cur_pd;
                double best_tx = sh.get_translate_double_x();
                double best_ty = sh.get_translate_double_y();

                // 计算形状的包围半径（用于圆内采样）
                double max_sq = 0.0;
                auto first_v = sh.base->outer_boundary().vertices_begin();
                for (auto v = sh.base->outer_boundary().vertices_begin(); v != sh.base->outer_boundary().vertices_end(); ++v) {
                    double dx = CGAL::to_double(v->x() - first_v->x());
                    double dy = CGAL::to_double(v->y() - first_v->y());
                    double d2 = dx*dx + dy*dy;
                    if (d2 > max_sq) max_sq = d2;
                }
                double rB = std::sqrt(max_sq);

                if (!layout.sheets.empty() && layout.sheets[0].is_circle()) {
                    geo::FT Rft = layout.sheets[0].get_radius();
                    double R = CGAL::to_double(Rft);
                    double innerR = R - rB;
                    if (innerR <= 0) {
                        // 无法放入，移动到远处并标记为未放置
                        sh.set_translate(-1e6, -1e6);
                        layout.mark_part_unplaced(i);
                        any_improved = true;
                        continue;
                    }

                    // 在圆内按面积均匀采样候选点
                    for (int t = 0; t < TRIES; ++t) {
                        double u = uniform01();
                        double rr = innerR * std::sqrt(u);
                        double theta = 2.0 * PI * uniform01();
                        double cx = R + rr * std::cos(theta);
                        double cy = R + rr * std::sin(theta);
                        double tx = cx - CGAL::to_double(first_v->x());
                        double ty = cy - CGAL::to_double(first_v->y());

                        sh.set_translate(tx, ty);
                        if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;

                        double total_pd = 0.0;
                        // 计算与已放置其他零件的pd之和（可剪枝）
                        for (size_t j = 0; j < layout.poly_num; ++j) {
                            if (i == j) continue;
                            auto& other = layout.sheet_parts[0][j];
                            double ox = other.get_translate_double_x();
                            double oy = other.get_translate_double_y();
                            if (ox < -1e5 || oy < -1e5) continue; // 未放置
                            auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                                sh.base, sh.get_rotation(), sh.allowed_rotations, layout);
                            FT pd = comp_pd(nfp, tx - ox, ty - oy, layout);
                            total_pd += CGAL::to_double(pd);
                            if (total_pd >= best_pd) break; // 剪枝
                        }
                        if (total_pd < best_pd) {
                            best_pd = total_pd; best_tx = tx; best_ty = ty;
                        }
                    }
                } else {
                    // 矩形板材随机采样
                    double height = CGAL::to_double(layout.sheets[0].get_height());
                    double length = CGAL::to_double(layout.cur_length);
                    if (length <= 0 || height <= 0) continue;
                    for (int t = 0; t < TRIES; ++t) {
                        double tx = uniform01() * length;
                        double ty = uniform01() * height;
                        sh.set_translate(tx, ty);

                        double total_pd = 0.0;
                        for (size_t j = 0; j < layout.poly_num; ++j) {
                            if (i == j) continue;
                            auto& other = layout.sheet_parts[0][j];
                            double ox = other.get_translate_double_x();
                            double oy = other.get_translate_double_y();
                            if (ox < -1e5 || oy < -1e5) continue;
                            auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                                sh.base, sh.get_rotation(), sh.allowed_rotations, layout);
                            FT pd = comp_pd(nfp, tx - ox, ty - oy, layout);
                            total_pd += CGAL::to_double(pd);
                            if (total_pd >= best_pd) break;
                        }
                        if (total_pd < best_pd) { best_pd = total_pd; best_tx = tx; best_ty = ty; }
                    }
                }

                // 接受改进
                if (best_pd < cur_pd) {
                    sh.set_translate(best_tx, best_ty);
                    // 更新pd表
                    for (size_t j = 0; j < layout.poly_num; ++j) {
                        if (i == j) continue;
                        auto& other = layout.sheet_parts[0][j];
                        double ox = other.get_translate_double_x();
                        double oy = other.get_translate_double_y();
                        auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                            sh.base, sh.get_rotation(), sh.allowed_rotations, layout);
                        FT pd = comp_pd(nfp, best_tx - ox, best_ty - oy, layout);
                        layout.set_pd(i, j, pd);
                    }
                    any_improved = true;
                }
            }

            double cur_total = CGAL::to_double(layout.get_pure_total_pd());
            if (cur_total < geo::BIAS) return true; // 已经无重叠
            if (!any_improved) break; // 局部无改进，则停止

            ++iterations;
        }

        return false;
    }

    // 尝试贪婪地放置一个零件（检查是否能放置）
    bool try_place_part_greedy(Layout& layout, const Polygon_with_holes_2& part, double /*R*/) {
        if (layout.sheets.empty() || !layout.sheets[0].is_circle()) return false;  // 检查是否为圆形板材

        geo::FT Rft = layout.sheets[0].get_radius();  // 获取板材半径
        double R = CGAL::to_double(Rft);
        const int RINGS = 6;  // 圆环数
        const int ANGULAR = 32;  // 角度采样数

        geo::Polygon_2 outer = part.outer_boundary();  // 获取零件外边界
        auto first_v = outer.vertices_begin();  // 第一个顶点

        // 计算零件的最大半径
        double max_sq = 0.0;
        for (auto v = outer.vertices_begin(); v != outer.vertices_end(); ++v) {
            double dx = CGAL::to_double(v->x() - first_v->x());
            double dy = CGAL::to_double(v->y() - first_v->y());
            double d2 = dx * dx + dy * dy;
            if (d2 > max_sq) max_sq = d2;
        }

        double rB = std::sqrt(max_sq);  // 零件的最大半径
        double innerR = R - rB;  // 内接圆半径
        if (innerR <= 0.0) return false;  // 如果内接圆半径非正，无法放置

        // 在多个圆环上采样候选位置
        for (int ring = 0; ring <= RINGS; ++ring) {
            double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);  // 当前圆环半径
            int samples = (ring == 0) ? 1 : ANGULAR;  // 角度采样数

            for (int s = 0; s < samples; ++s) {
                double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);  // 角度
                double centre_x = R + ring_r * std::cos(theta);  // 候选位置中心x坐标
                double centre_y = R + ring_r * std::sin(theta);  // 候选位置中心y坐标

                bool inside = true;  // 标记是否所有顶点都在圆内

                // 检查零件的所有顶点是否都在圆内
                for (auto v = outer.vertices_begin(); v != outer.vertices_end(); ++v) {
                    double dx = CGAL::to_double(v->x()) - centre_x;  // 顶点到圆心的x距离
                    double dy = CGAL::to_double(v->y()) - centre_y;  // 顶点到圆心的y距离
                    // 如果距离平方大于半径平方（加上容差），则不在圆内
                    if (dx * dx + dy * dy > R * R + 1e-9) {
                        inside = false;
                        break;
                    }
                }

                // 如果所有顶点都在圆内，返回成功
                if (inside) return true;
            }
        }

        return false;  // 所有候选位置都不满足条件，返回失败
    }

    // 计算当前利用率
    static double compute_current_utilization(const Layout& layout) {
        if (layout.sheets.empty() || !layout.sheets[0].is_circle()) return 0.0;  // 检查是否为圆形板材

        geo::FT R = layout.sheets[0].get_radius();  // 获取板材半径
        double parts_area = 0.0;  // 零件总面积

        // 遍历所有零件，累加完全在圆内的零件面积
        for (size_t i = 0; i < layout.poly_num; ++i) {
            const auto& p = layout.sheet_parts[0][i];  // 获取零件
            if (all_vertices_inside_circle(p.transformed, R))  // 如果零件完全在圆内
                parts_area += CGAL::to_double(geo::pwh_area(p.transformed));  // 累加面积
        }

        double Rd = CGAL::to_double(R);  // 半径转换为double
        double sheet_area = PI * Rd * Rd;  // 计算圆形板材面积
        if (sheet_area <= 0.0) return 0.0;  // 面积无效则返回0

        double util = parts_area / sheet_area;  // 计算利用率
        if (util > 1.0) util = 1.0;  // 限制利用率不超过1.0
        return util;  // 返回利用率
    }

    // 随机扰动一个零件的位置
    static void perturb_random_part(Layout& layout) {
        if (layout.sheet_parts.empty() || layout.sheet_parts[0].empty()) return;  // 没有零件则返回

        size_t n = layout.sheet_parts[0].size();  // 零件数量
        size_t idx = (size_t)(std::rand() % n);  // 随机选择一个零件
        auto& sh = layout.sheet_parts[0][idx];  // 获取零件引用

        // 生成随机扰动
        double jitter = 1.0 + 10.0 * ((double)std::rand() / (double)RAND_MAX);  // 扰动幅度
        double ang = 2.0 * PI * ((double)std::rand() / (double)RAND_MAX);  // 扰动方向
        double tx = sh.get_translate_double_x() + jitter * std::cos(ang);  // 新x坐标
        double ty = sh.get_translate_double_y() + jitter * std::sin(ang);  // 新y坐标

        sh.set_translate(tx, ty);  // 更新零件位置
    }

    // 尝试添加新零件（委托给attempt_place_unplaced_part）
    static void try_add_new_part(Layout& layout) { attempt_place_unplaced_part(layout); }

    // 尝试移除最小的零件
    static void try_remove_small_part(Layout& layout) {
        if (layout.sheet_parts.empty() || layout.sheet_parts[0].empty()) return;  // 没有零件则返回

        size_t n = layout.sheet_parts[0].size();  // 零件数量
        size_t min_idx = 0;  // 最小零件的索引
        double min_a = std::numeric_limits<double>::infinity();  // 最小零件面积

        // 查找面积最小的零件
        for (size_t i = 0; i < n; ++i) {
            double a = CGAL::to_double(geo::pwh_area(layout.sheet_parts[0][i].transformed));  // 零件面积
            if (a < min_a) {
                min_a = a;  // 更新最小面积
                min_idx = i;  // 更新最小零件索引
            }
        }

        const double FAR = -1e6;  // 远距离坐标
        // 将最小零件移到远处
        layout.sheet_parts[0][min_idx].set_translate(FAR, FAR);
        layout.mark_part_unplaced(min_idx);  // 标记为未放置
    }

    // 模拟退火优化
    static void simulated_annealing_optimization(Layout& layout,
        size_t max_time,  // 最大运行时间
        const clock_t& start,  // 开始时间
        std::function<void(const Solution&)> ProgressHandler,  // 进度回调函数
        volatile bool* requestQuit) {  // 退出标志

        double initial_temperature = 1.0;
        double final_temperature = 0.01;
        double cooling_rate = 0.95;

        std::mt19937_64 rng((unsigned long)clock() ^ (unsigned long)std::hash<std::thread::id>()(std::this_thread::get_id()));
        auto random_double = [&](double lo = 0.0, double hi = 1.0) {
            std::uniform_real_distribution<double> dist(lo, hi);
            return dist(rng);
        };

        Layout current_layout = layout;
        Layout best_layout = layout;
        double best_util = compute_current_utilization(layout);
        double current_util = best_util;
        double temperature = initial_temperature;

        // 局部优化：尝试优化单个零件以减少其穿透量
        auto optimize_single_part = [&](Layout& L, size_t part_idx) -> bool {
            if (L.sheet_parts.empty() || part_idx >= L.poly_num) return false;
            auto& sh = L.sheet_parts[0][part_idx];
            double cur_pd = CGAL::to_double(L.get_one_polygon_pd(part_idx));
            if (cur_pd <= geo::BIAS) return false;

            const int TRIES = 10;
            double best_pd = cur_pd;
            double best_tx = sh.get_translate_double_x();
            double best_ty = sh.get_translate_double_y();

            // 计算包围半径
            auto first_v = sh.base->outer_boundary().vertices_begin();
            double max_sq = 0.0;
            for (auto v = sh.base->outer_boundary().vertices_begin(); v != sh.base->outer_boundary().vertices_end(); ++v) {
                double dx = CGAL::to_double(v->x() - first_v->x());
                double dy = CGAL::to_double(v->y() - first_v->y());
                double d2 = dx*dx + dy*dy;
                if (d2 > max_sq) max_sq = d2;
            }
            double rB = std::sqrt(max_sq);

            if (!L.sheets.empty() && L.sheets[0].is_circle()) {
                geo::FT Rft = L.sheets[0].get_radius();
                double R = CGAL::to_double(Rft);
                double innerR = R - rB;
                if (innerR <= 0) return false;

                for (int t = 0; t < TRIES; ++t) {
                    double u = random_double();
                    double rr = innerR * std::sqrt(u);
                    double theta = 2.0 * PI * random_double();
                    double cx = R + rr * std::cos(theta);
                    double cy = R + rr * std::sin(theta);
                    double tx = cx - CGAL::to_double(first_v->x());
                    double ty = cy - CGAL::to_double(first_v->y());

                    sh.set_translate(tx, ty);
                    if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;

                    double total_pd = 0.0;
                    for (size_t j = 0; j < L.poly_num; ++j) {
                        if (j == part_idx) continue;
                        auto& other = L.sheet_parts[0][j];
                        double ox = other.get_translate_double_x();
                        double oy = other.get_translate_double_y();
                        if (ox < -1e5 || oy < -1e5) continue;

                        auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                            sh.base, sh.get_rotation(), sh.allowed_rotations, L);
                        FT pd = comp_pd(nfp, tx - ox, ty - oy, L);
                        total_pd += CGAL::to_double(pd);
                        if (total_pd >= best_pd) break;
                    }

                    if (total_pd < best_pd) { best_pd = total_pd; best_tx = tx; best_ty = ty; }
                }
            } else {
                // 矩形板材小幅扰动
                double base_x = sh.get_translate_double_x();
                double base_y = sh.get_translate_double_y();
                double radius = 5.0;
                for (int t = 0; t < TRIES; ++t) {
                    double ang = 2.0 * PI * random_double();
                    double r = radius * random_double();
                    double tx = base_x + r * std::cos(ang);
                    double ty = base_y + r * std::sin(ang);
                    sh.set_translate(tx, ty);

                    double total_pd = 0.0;
                    for (size_t j = 0; j < L.poly_num; ++j) {
                        if (j == part_idx) continue;
                        auto& other = L.sheet_parts[0][j];
                        double ox = other.get_translate_double_x();
                        double oy = other.get_translate_double_y();
                        if (ox < -1e5 || oy < -1e5) continue;
                        auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                            sh.base, sh.get_rotation(), sh.allowed_rotations, L);
                        FT pd = comp_pd(nfp, tx - ox, ty - oy, L);
                        total_pd += CGAL::to_double(pd);
                        if (total_pd >= best_pd) break;
                    }
                    if (total_pd < best_pd) { best_pd = total_pd; best_tx = tx; best_ty = ty; }
                }
            }

            if (best_pd < cur_pd) { sh.set_translate(best_tx, best_ty); return true; }
            return false;
        };

        // 交换两个零件的位置
        auto swap_two_parts = [&](Layout& L) {
            if (L.sheet_parts.empty() || L.sheet_parts[0].size() < 2) return;
            size_t n = L.sheet_parts[0].size();
            size_t a = (size_t)std::min<double>(std::floor(random_double() * n), (double)n - 1);
            size_t b = (size_t)std::min<double>(std::floor(random_double() * n), (double)n - 1);
            if (a == b) return;
            auto& A = L.sheet_parts[0][a]; auto& B = L.sheet_parts[0][b];
            double ax = A.get_translate_double_x(); double ay = A.get_translate_double_y();
            double bx = B.get_translate_double_x(); double by = B.get_translate_double_y();
            A.set_translate(bx, by); B.set_translate(ax, ay);
        };

        // 主模拟退火循环
        while (((double)(clock() - start) / CLOCKS_PER_SEC) < max_time) {
            if (*requestQuit) return;

            for (int i = 0; i < 100; ++i) {
                if (*requestQuit) return;
                Layout new_layout = current_layout;
                bool changed = false;

                int operation = (int)std::floor(random_double() * 6.0);
                switch (operation) {
                    case 0: perturb_random_part(new_layout); changed = true; break;
                    case 1: swap_two_parts(new_layout); changed = true; break;
                    case 2: try_add_new_part(new_layout); changed = true; break;
                    case 3: try_remove_small_part(new_layout); changed = true; break;
                    case 4: attempt_place_unplaced_part(new_layout); changed = true; break;
                    case 5: {
                        if (new_layout.poly_num > 0) {
                            size_t idx = (size_t)std::min<double>(std::floor(random_double() * new_layout.poly_num), (double)new_layout.poly_num - 1);
                            changed = optimize_single_part(new_layout, idx);
                        }
                        break;
                    }
                    default: perturb_random_part(new_layout); changed = true; break;
                }

                if (!changed) continue;

                double new_util = compute_current_utilization(new_layout);
                double delta = current_util - new_util;

                bool accept = false;
                if (delta < 0.0) accept = true;
                else if (temperature > 0.0) {
                    double prob = std::exp(-delta / temperature);
                    if (prob > random_double()) accept = true;
                }

                if (accept) {
                    current_layout = new_layout; current_util = new_util;
                    if (new_util > best_util) {
                        best_util = new_util; best_layout = new_layout;
                        double disp_len = 0.0; if (!best_layout.sheets.empty()) disp_len = CGAL::to_double(best_layout.sheets[0].get_diameter());
                        ProgressHandler(Solution(disp_len, best_util, ((double)(clock() - start) / CLOCKS_PER_SEC), best_layout.sheet_parts[0]));
                    }
                }
            }

            temperature *= cooling_rate;
            if (temperature < final_temperature) temperature = initial_temperature * 0.5;
        }

        layout = best_layout;
        layout.best_utilization = best_util;
    }

    // GOMH（Guided Online Meta-Heuristic）算法主函数
    void GOMH(Layout& layout,
        size_t max_time,  // 最大运行时间
        std::function<void(const Solution&)> ProgressHandler,  // 进度回调函数
        volatile bool* requestQuit) {  // 退出标志

        clock_t start = clock();
        double fixed_radius = 150.0;

        // 设置圆形板材半径
        if (!layout.sheets.empty()) layout.sheets[0].set_radius(geo::FT(fixed_radius));

        // 获取贪婪初始解
        get_greedy_initial_solution(layout);

        // 计算显示长度
        double length_display = CGAL::to_double(layout.best_length);
        if (!layout.sheets.empty() && layout.sheets[0].is_circle()) {
            length_display = CGAL::to_double(layout.sheets[0].get_diameter());
        }

        double current_util = compute_current_utilization(layout);
        layout.best_utilization = current_util;
        layout.best_result = layout.sheet_parts;

        // 报告初始解
        ProgressHandler(Solution(length_display, layout.best_utilization,
            ((double)(clock() - start) / CLOCKS_PER_SEC), layout.best_result[0]));

        // 收缩布局（可选）
        shrink(layout);

        // 直接调用模拟退火优化，不再使用主循环
        simulated_annealing_optimization(layout, max_time, start, ProgressHandler, requestQuit);

        // 最终收缩并更新最佳利用率
        shrink(layout);
        layout.best_utilization = compute_current_utilization(layout);
    }

    // 贪婪获取初始解
    static void get_greedy_initial_solution(Layout& layout) {
        if (layout.sheets.empty() || !layout.sheets[0].is_circle())
            throw std::runtime_error("get_greedy_initial_solution requires a circular sheet");  // 检查是否为圆形板材

        geo::FT Rft = layout.sheets[0].get_radius();  // 获取板材半径
        double R = CGAL::to_double(Rft);
        double center_x = R, center_y = R;  // 圆心坐标

        size_t total_parts = layout.sheet_parts.size() ? layout.sheet_parts[0].size() : 0;  // 零件总数
        if (total_parts == 0) return;  // 没有零件则返回

        // 存储零件索引和面积的结构体
        struct PartArea { size_t idx; double area; };
        std::vector<PartArea> parts;
        parts.reserve(total_parts);

        // 计算所有零件的面积并存储
        for (size_t i = 0; i < total_parts; ++i)
            parts.push_back({ i, CGAL::to_double(geo::pwh_area(layout.sheet_parts[0][i].transformed)) });

        // 按面积从大到小排序（大零件优先放置）
        std::sort(parts.begin(), parts.end(), [](const PartArea& A, const PartArea& B) { return A.area > B.area; });

        std::vector<char> placed(total_parts, 0);  // 标记零件是否已放置（0未放置，1已放置）
        const int RINGS = 8;  // 圆环数
        const int ANGULAR = 64;  // 角度采样数

        // 按面积从大到小遍历零件
        for (const auto& pa : parts) {
            size_t i = pa.idx;  // 零件索引
            auto& sh = layout.sheet_parts[0][i];  // 获取零件引用
            uint32_t rotation_i = sh.get_rotation();  // 获取零件旋转角度
            auto rotate = geo::get_rotate(rotation_i, sh.allowed_rotations);  // 获取旋转变换

            // 应用旋转变换
            Polygon_with_holes_2 rotated = rotate ? geo::transform_polygon_with_holes(*rotate, *sh.base) : *sh.base;
            auto first_v = rotated.outer_boundary().vertices_begin();  // 第一个顶点

            // 计算零件的最大半径
            double max_sq = 0.0;
            for (auto vit = rotated.outer_boundary().vertices_begin(); vit != rotated.outer_boundary().vertices_end(); ++vit) {
                double dx = CGAL::to_double(vit->x() - first_v->x());
                double dy = CGAL::to_double(vit->y() - first_v->y());
                double d2 = dx * dx + dy * dy;
                if (d2 > max_sq) max_sq = d2;
            }

            double rB = std::sqrt(max_sq);  // 零件的最大半径
            double innerR = R - rB;  // 内接圆半径
            if (innerR <= 0.0) continue;  // 如果内接圆半径非正，跳过该零件

            bool success = false;  // 标记是否成功放置

            // 在多个圆环上采样候选位置
            for (int ring = 0; ring <= RINGS && !success; ++ring) {
                double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);  // 当前圆环半径
                int samples = (ring == 0) ? 1 : ANGULAR;  // 角度采样数

                for (int s = 0; s < samples; ++s) {
                    double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);  // 角度
                    double cx = center_x + ring_r * std::cos(theta);  // 候选位置中心x坐标
                    double cy = center_y + ring_r * std::sin(theta);  // 候选位置中心y坐标

                    // 计算平移量
                    double tx = cx - CGAL::to_double(first_v->x());
                    double ty = cy - CGAL::to_double(first_v->y());
                    sh.set_translate(tx, ty);  // 设置零件位置

                    // 检查零件是否完全在圆内
                    if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;

                    double overlap_sum = 0.0;  // 总重叠量

                    // 检查与所有已放置零件的重叠情况
                    for (size_t j = 0; j < total_parts; ++j) {
                        if (!placed[j]) continue;  // 跳过未放置的零件

                        auto& other = layout.sheet_parts[0][j];  // 已放置的零件
                        double ox = other.get_translate_double_x();  // 已放置零件的x坐标
                        double oy = other.get_translate_double_y();  // 已放置零件的y坐标

                        // 计算NFP
                        auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                            sh.base, rotation_i, sh.allowed_rotations, layout);
                        // 计算穿透距离
                        FT pd = comp_pd(nfp, tx - ox, ty - oy, layout);
                        overlap_sum += CGAL::to_double(pd);  // 累加穿透距离

                        // 如果有重叠，提前退出
                        if (overlap_sum > geo::BIAS) break;
                    }

                    // 如果没有重叠，标记为已放置
                    if (overlap_sum <= geo::BIAS) {
                        placed[i] = 1;
                        success = true;
                        break;
                    }
                }
            }
        }

        const double FAR = -1e6;  // 远距离坐标

        // 将未放置的零件移到远处
        for (size_t i = 0; i < total_parts; ++i) {
            if (!placed[i]) layout.sheet_parts[0][i].set_translate(FAR, FAR);
        }

        // 计算所有已放置零件之间的穿透距离
        for (size_t i = 0; i < total_parts; ++i) {
            if (!placed[i]) continue;  // 跳过未放置的零件

            for (size_t j = i + 1; j < total_parts; ++j) {
                if (!placed[j]) continue;  // 跳过未放置的零件

                auto& A = layout.sheet_parts[0][i];  // 零件A
                auto& B = layout.sheet_parts[0][j];  // 零件B

                double ax = A.get_translate_double_x();  // A的x坐标
                double ay = A.get_translate_double_y();  // A的y坐标
                double bx = B.get_translate_double_x();  // B的x坐标
                double by = B.get_translate_double_y();  // B的y坐标

                // 计算NFP
                auto& nfp = comp_nfp(A.base, A.get_rotation(), A.allowed_rotations,
                    B.base, B.get_rotation(), B.allowed_rotations, layout);
                // 计算穿透距离
                FT pd = comp_pd(nfp, ax - bx, ay - by, layout);
                layout.set_pd(i, j, pd);  // 设置穿透距离
            }
        }

        // 统计信息
        size_t count = 0;  // 已放置零件计数
        double totalArea = 0.0;  // 已放置零件总面积

        for (size_t i = 0; i < total_parts; ++i) {
            if (placed[i]) {
                ++count;
                totalArea += CGAL::to_double(geo::pwh_area(layout.sheet_parts[0][i].transformed));
            }
        }

        // 输出贪婪算法的结果
        std::clog << "Greedy placed " << count << " of " << total_parts << " parts, total area " << totalArea << std::endl;

        layout.update_cur_length();  // 更新当前长度
        layout.best_length = layout.cur_length;  // 设置最佳长度
        layout.best_result = layout.sheet_parts;  // 设置最佳结果
        layout.best_utilization = compute_util(layout, 0.0);  // 计算最佳利用率
    }

} // namespace nesting