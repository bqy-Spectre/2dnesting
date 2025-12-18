#pragma once
#include "algorithm.h"
#include "hasher.h"
#include "lru_size.h"
#include "part.h"
#include "sheet.h"
#include "transformed_shape.h"
namespace nesting {

    typedef typename geo::FT FT;
    typedef typename geo::Polygon_with_holes_2 Polygon_with_holes_2;
    typedef typename geo::Transformation Transformation;
    typedef typename geo::Vector_2 Vector_2;
    typedef typename hash::NFPCacheKey NFPCacheKey;
    typedef typename hash::NFPCacheValue NFPCacheValue;
    typedef typename hash::NFPCacheKeyHasher NFPCacheKeyHasher;

    // 算法运行时的布局情况，用于缓存已计算的NFP和已放置的多边形
    struct Layout {
        Layout(std::vector<Item> items, std::vector<Sheet> sheets, FT _area = -1);

    private:
        // 当前布局下的pd，压缩到一维。pd[i][j]映射到pd[n*i-i*(i+1)/2+(j-i-1)]
        std::vector<FT> glob_pd;
        std::vector<double> miu;

    public:
        // 超参数
        double rinc{ 0.01 };
        double rdec{ 0.04 };
        size_t maxIterations{ 50 };
        uint32_t max_pd_cache{ 200000 };
        // 多边形总数
        size_t poly_num{ 0 };
        size_t pd_count{ 0 };
        size_t pd_miss{ 0 };
        // 多板，板可能为矩形或圆形
        std::vector<Sheet> sheets;
        // 已计算的nfp缓存
        std::unordered_map<NFPCacheKey, NFPCacheValue, NFPCacheKeyHasher> nfp_cache;
        // 已计算的规范多边形集合
        std::unordered_map<Polygon_with_holes_2,
            std::shared_ptr<Polygon_with_holes_2>,
            hash::PolygonHasher>
            canonical_polygons;
        // LRU PD缓存，线程不安全
        emlru_size::lru_cache<hash::PDCacheKey, FT, hash::PDCacheKeyHasher> pd_cache;
        // sheet_parts[i]代表sheet i上存储的变形shape
        std::vector<std::vector<TransformedShape>> sheet_parts;
        std::vector<std::vector<TransformedShape>> best_result;
        // 当前长度
        FT cur_length{ -1 };
        // 最短长度
        FT best_length{ -1 };
        // 最佳利用率
        double best_utilization{ 0 };
        // 排样长度的下界
        FT lower_length{ -1 };
        // 总面积
        FT area{ 0 };
        bool isRequestQuit{ false };

        // --- 新增：零件管理字段与方法 ---
        struct ManagedPart {
            Polygon_with_holes_2 polygon;
            double area{0.0};
            bool placed{false};
            size_t original_index{0};
        };

        std::vector<ManagedPart> all_parts;   // 所有可用零件
        std::vector<size_t> placed_indices;   // 已放置零件的索引（索引指向 all_parts）
        double total_parts_area{0.0};         // 已放置零件的总面积

        // 初始化 all_parts 列表
        void init_all_parts(const std::vector<Polygon_with_holes_2>& parts) {
            all_parts.clear();
            all_parts.reserve(parts.size());
            for (size_t i = 0; i < parts.size(); ++i) {
                ManagedPart mp;
                mp.polygon = parts[i];
                mp.area = CGAL::to_double(geo::pwh_area(parts[i]));
                mp.placed = false;
                mp.original_index = i;
                all_parts.push_back(std::move(mp));
            }
            placed_indices.clear();
            total_parts_area = 0.0;
        }

        // 返回未放置零件的指针列表
        std::vector<ManagedPart*> get_unplaced_parts() {
            std::vector<ManagedPart*> result;
            result.reserve(all_parts.size());
            for (auto& mp : all_parts) {
                if (!mp.placed) result.push_back(&mp);
            }
            return result;
        }

        // 标记某个 ManagedPart 为已放置（并更新 placed_indices 与 total_parts_area）
        void mark_part_placed(size_t all_parts_index, size_t placed_index_in_sheet_parts = SIZE_MAX) {
            if (all_parts_index >= all_parts.size()) return;
            if (all_parts[all_parts_index].placed) return;
            all_parts[all_parts_index].placed = true;
            placed_indices.push_back(all_parts_index);
            total_parts_area += all_parts[all_parts_index].area;
            // placed_index_in_sheet_parts 可用于将 all_parts 索引映射到 sheet_parts 的位置（如需要）
        }

        // 标记某个 ManagedPart 为未放置（并更新 placed_indices 与 total_parts_area）
        void mark_part_unplaced(size_t all_parts_index) {
            if (all_parts_index >= all_parts.size()) return;
            if (!all_parts[all_parts_index].placed) return;
            all_parts[all_parts_index].placed = false;
            // 从 placed_indices 中移除该索引
            auto it = std::find(placed_indices.begin(), placed_indices.end(), all_parts_index);
            if (it != placed_indices.end()) placed_indices.erase(it);
            total_parts_area -= all_parts[all_parts_index].area;
            if (total_parts_area < 1e-12) total_parts_area = 0.0;
        }

        // --- 其它现有声明保持不变 ---
        /// @brief 存放一个规范多边形，规范多边形的意思是其第一个顶点平移到(0, 0)处
        /// @param poly 需要存放的规范多边形
        /// @return 指向该多边形的一个指针
        Polygon_with_holes_2* get_canonical_polygon(const Polygon_with_holes_2& poly);

        // 设置polygon[i]和polygon[j]之间的pd
        inline void set_pd(size_t i, size_t j, FT pd) {
            if (i == j) {
                return;
            }
            if (i > j) {
                i ^= j;
                j ^= i;
                i ^= j;
            }
            size_t idx = poly_num * i - i * (i + 1) / 2 + (j - i - 1);
            glob_pd[idx] = pd;
        }

        // 获取polygon[i]和polygon[j]之间的加权pd
        inline FT get_pd(size_t i, size_t j) {
            if (i == j) {
                return 0;
            }
            if (i > j) {
                i ^= j;
                j ^= i;
                i ^= j;
            }
            size_t idx = poly_num * i - i * (i + 1) / 2 + (j - i - 1);
            return glob_pd[idx] * miu[idx];
        }

        // 初始化权重均为1
        inline void initialize_miu() {
            auto size = miu.size();
            for (size_t i = 0; i < size; i++) {
                miu[i] = 1;
            }
        }

        // 更新权重
        void update_miu();

        // 获取权重
        inline double get_miu(size_t i, size_t j) {
            if (i == j) {
                return 0;
            }
            else if (i > j) {
                i ^= j;
                j ^= i;
                i ^= j;
            }
            size_t idx = poly_num * i - i * (i + 1) / 2 + (j - i - 1);
            return miu[idx];
        }

        // 获取关于polygon[p]的加权pd之和
        inline FT get_one_polygon_pd(size_t p) {
            FT t = 0;
            for (size_t i = 0; i < poly_num; i++) {
                t += get_pd(p, i);
            }
            return t;
        }

        // 获取所有polygon的加权pd之和
        inline FT get_total_pd() {
            FT total_pd = 0;
            size_t size = glob_pd.size();
            for (size_t i = 0; i < size; i++) {
                total_pd += glob_pd[i] * miu[i];
            }
            return total_pd;
        }

        // 获取所有polygon的纯pd之和
        inline FT get_pure_total_pd() {
            FT total_pd = 0;
            size_t size = glob_pd.size();
            for (size_t i = 0; i < size; i++) {
                total_pd += glob_pd[i];
            }
            return total_pd;
        }

        // 更新现有长度
        void update_cur_length();

        void debug_glob_pd();
    };
}  // namespace nesting