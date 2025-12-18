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

// Define PI for MSVC where M_PI may be unavailable
#ifndef PI
static constexpr double PI = 3.14159265358979323846;
#endif

namespace nesting {

using Polygon_2 = geo::Polygon_2;

// Forward declarations
static bool all_vertices_inside_circle(const Polygon_with_holes_2& pwh, const geo::FT& R);
static void get_greedy_initial_solution(Layout& layout);
static void attempt_place_unplaced_part(Layout& L);
bool minimize_overlap(Layout& layout, volatile bool* requestQuit);

// Compute utilization: circular uses parts fully inside circle; rectangular uses placed parts area
static double compute_util(const Layout& layout, double /*length*/) {
    if (layout.sheets.empty()) return 0.0;
    if (layout.sheets[0].is_circle()) {
        geo::FT R = layout.sheets[0].get_radius();
        if (R <= geo::FT(0)) return 0.0;
        double parts_area = 0.0;
        if (!layout.sheet_parts.empty()) {
            for (size_t i = 0; i < layout.sheet_parts[0].size(); ++i) {
                const auto& p = layout.sheet_parts[0][i];
                if (all_vertices_inside_circle(p.transformed, R)) {
                    parts_area += CGAL::to_double(geo::pwh_area(p.transformed));
                }
            }
        }
        double Rd = CGAL::to_double(R);
        double sheet_area = PI * Rd * Rd;
        if (sheet_area <= 0.0) return 0.0;
        double util = parts_area / sheet_area;
        if (util > 1.0) util = 1.0;
        return util;
    }
    // rectangular fallback
    double parts_area = 0.0;
    const double FAR_THRESHOLD = -1e5;
    if (!layout.sheet_parts.empty()) {
        for (size_t i = 0; i < layout.sheet_parts[0].size(); ++i) {
            const auto& p = layout.sheet_parts[0][i];
            double tx = p.get_translate_double_x();
            double ty = p.get_translate_double_y();
            if (tx > FAR_THRESHOLD && ty > FAR_THRESHOLD) {
                parts_area += CGAL::to_double(geo::pwh_area(p.transformed));
            }
        }
    }
    double height = CGAL::to_double(layout.sheets[0].get_height());
    double length = CGAL::to_double(layout.cur_length);
    if (length <= 0.0 || height <= 0.0) return 0.0;
    const double EPS = 1e-12;
    double rect_area = length * height + EPS;
    double util = parts_area / rect_area;
    if (util > 1.0) util = 1.0;
    return util;
}

static bool all_vertices_inside_circle(const Polygon_with_holes_2& pwh, const geo::FT& R) {
    if (R <= geo::FT(0)) return false;
    double Rd = CGAL::to_double(R);
    double R2 = Rd * Rd;
    for (auto vit = pwh.outer_boundary().vertices_begin(); vit != pwh.outer_boundary().vertices_end(); ++vit) {
        double dx = CGAL::to_double(vit->x()) - Rd;
        double dy = CGAL::to_double(vit->y()) - Rd;
        if (dx*dx + dy*dy > R2 + 1e-12) return false;
    }
    for (auto hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit) {
        for (auto vit = hit->vertices_begin(); vit != hit->vertices_end(); ++vit) {
            double dx = CGAL::to_double(vit->x()) - Rd;
            double dy = CGAL::to_double(vit->y()) - Rd;
            if (dx*dx + dy*dy > R2 + 1e-12) return false;
        }
    }
    return true;
}

NFPCacheValue& comp_nfp(const Polygon_with_holes_2* poly_A,
    const uint32_t rotation_A,
    const uint32_t allowed_rotation_A,
    const Polygon_with_holes_2* poly_B,
    const uint32_t rotation_B,
    const uint32_t allowed_rotation_B,
    Layout& layout) {
    NFPCacheKey key(poly_A, poly_B, rotation_A, rotation_B);
    auto kv = layout.nfp_cache.find(key);
    if (kv == layout.nfp_cache.end()) {
        Transformation scale(CGAL::SCALING, -1);
        auto rotate_A = geo::get_rotate(rotation_A, allowed_rotation_A);
        auto rotate_B = geo::get_rotate(rotation_B, allowed_rotation_B);
        Polygon_with_holes_2 minus_B;
        Polygon_with_holes_2 rotated_A;
        if (rotate_B) minus_B = geo::transform_polygon_with_holes(scale * (*rotate_B), *poly_B);
        else minus_B = geo::transform_polygon_with_holes(scale, *poly_B);
        if (rotate_A) rotated_A = geo::transform_polygon_with_holes(*rotate_A, *poly_A);
        else rotated_A = *poly_A;
        Polygon_with_holes_2 nfp = CGAL::minkowski_sum_2(rotated_A, minus_B);
        geo::strict_simplify(nfp);
        auto bbox = nfp.bbox();
        NFPCacheValue v;
        v.nfp = nfp;
        v.xmin = CGAL::to_double(bbox.xmin());
        v.xmax = CGAL::to_double(bbox.xmax());
        v.ymin = CGAL::to_double(bbox.ymin());
        v.ymax = CGAL::to_double(bbox.ymax());
        auto _kv = layout.nfp_cache.emplace(key, v);
        return _kv.first->second;
    }
    return kv->second;
}

FT comp_pd(NFPCacheValue& v, double px, double py, Layout& layout) {
    auto& original_nfp = v.nfp;
    if (px <= v.xmin || px >= v.xmax || py <= v.ymin || py >= v.ymax) return geo::FT_ZERO;
    hash::PDCacheKey key(&original_nfp, px, py);
    FT pd;
    layout.pd_count++;
    auto iter = layout.pd_cache.find(key);
    if (iter != layout.pd_cache.end()) pd = iter->second;
    else {
        layout.pd_miss++;
        Point_2 relative_point(px, py);
        pd = geo::comp_pd(original_nfp, relative_point);
        layout.pd_cache.emplace(key, pd);
    }
    return pd;
}

void shrink(Layout& layout) {
    if (layout.sheets.empty()) return;
    if (!layout.sheets[0].is_circle()) return;
    std::vector<size_t> moved_indices;
    geo::FT R = layout.sheets[0].get_radius();
    double sheetR_d = CGAL::to_double(R);
    for (size_t i = 0; i < layout.poly_num; i++) {
        auto& p = layout.sheet_parts[0][i];
        if (!all_vertices_inside_circle(p.transformed, R)) {
            geo::Polygon_2 ifr_poly;
            auto first_v = p.transformed.outer_boundary().vertices_begin();
            double max_sq_d = 0.0;
            for (auto vtx = p.transformed.outer_boundary().vertices_begin(); vtx != p.transformed.outer_boundary().vertices_end(); ++vtx) {
                double dx = CGAL::to_double(vtx->x() - first_v->x());
                double dy = CGAL::to_double(vtx->y() - first_v->y());
                double d2 = dx * dx + dy * dy;
                if (d2 > max_sq_d) max_sq_d = d2;
            }
            geo::FT rB = geo::FT(std::sqrt(max_sq_d));
            geo::FT innerR = R - rB;
            if (innerR > geo::FT(0)) {
                const int SEGMENTS = 128;
                double cx = sheetR_d;
                double cy = sheetR_d;
                double rr = CGAL::to_double(innerR);
                for (int s = 0; s < SEGMENTS; ++s) {
                    double theta = 2.0 * PI * s / SEGMENTS;
                    double xpt = cx + rr * std::cos(theta);
                    double ypt = cy + rr * std::sin(theta);
                    ifr_poly.push_back(Point_2(xpt, ypt));
                }
                // translate polygon manually
                geo::Polygon_2 translated;
                for (auto vit = ifr_poly.vertices_begin(); vit != ifr_poly.vertices_end(); ++vit) {
                    double nx = CGAL::to_double(vit->x()) - CGAL::to_double(first_v->x());
                    double ny = CGAL::to_double(vit->y()) - CGAL::to_double(first_v->y());
                    translated.push_back(Point_2(nx, ny));
                }
                ifr_poly = translated;
            }
            if (ifr_poly.size() > 0) {
                auto bbox = ifr_poly.bbox();
                double ux = (double)std::rand() / (double)RAND_MAX;
                double uy = (double)std::rand() / (double)RAND_MAX;
                double random_x = bbox.x_span() * ux + bbox.xmin();
                double random_y = bbox.y_span() * uy + bbox.ymin();
                p.set_translate(random_x, random_y);
            } else {
                p.set_translate(CGAL::to_double(R), CGAL::to_double(R));
            }
            moved_indices.push_back(i);
        }
    }
    for (auto& i : moved_indices) {
        auto& p = layout.sheet_parts[0][i];
        auto double_px = CGAL::to_double(p.get_translate_ft_x());
        auto double_py = CGAL::to_double(p.get_translate_ft_y());
        for (size_t j = 0; j < layout.poly_num; j++) {
            if (i == j) continue;
            auto& q = layout.sheet_parts[0][j];
            auto double_qx = q.get_translate_double_x();
            auto double_qy = q.get_translate_double_y();
            auto& nfp = comp_nfp(q.base, q.get_rotation(), q.allowed_rotations, p.base, p.get_rotation(), p.allowed_rotations, layout);
            auto pd = comp_pd(nfp, double_px - double_qx, double_py - double_qy, layout);
            layout.set_pd(i, j, pd);
        }
    }
    layout.best_utilization = compute_util(layout, 0.0);
}

void get_init_solu(Layout& layout) {
    if (layout.sheets.empty() || !layout.sheets[0].is_circle()) throw std::runtime_error("get_init_solu requires a circular sheet");
    get_greedy_initial_solution(layout);
    if (layout.all_parts.empty()) {
        std::vector<Polygon_with_holes_2> parts;
        size_t total_parts = layout.sheet_parts.size() ? layout.sheet_parts[0].size() : 0;
        parts.reserve(total_parts);
        for (size_t i = 0; i < total_parts; ++i) parts.push_back(layout.sheet_parts[0][i].transformed);
        layout.init_all_parts(parts);
    }
    const double FAR = -1e6;
    layout.placed_indices.clear();
    layout.total_parts_area = 0.0;
    for (size_t i = 0; i < layout.poly_num && i < layout.all_parts.size(); ++i) {
        auto& sh = layout.sheet_parts[0][i];
        double tx = sh.get_translate_double_x();
        double ty = sh.get_translate_double_y();
        if (tx <= FAR/2 || ty <= FAR/2) layout.mark_part_unplaced(i);
        else {
            geo::FT R = layout.sheets[0].get_radius();
            if (all_vertices_inside_circle(sh.transformed, R)) layout.mark_part_placed(i);
            else layout.mark_part_unplaced(i);
        }
    }
    layout.update_cur_length();
    layout.best_length = layout.cur_length;
    layout.best_result = layout.sheet_parts;
    layout.best_utilization = compute_util(layout, 0.0);
}

static void attempt_place_unplaced_part(Layout& L) {
    if (L.sheets.empty() || !L.sheets[0].is_circle()) return;
    if (L.all_parts.empty()) return;
    auto unplaced = L.get_unplaced_parts();
    if (unplaced.empty()) return;
    size_t pick = (size_t)(std::rand() % unplaced.size());
    auto* mp = unplaced[pick];
    size_t idx = (size_t)(mp - &L.all_parts[0]);
    if (idx >= L.sheet_parts[0].size()) return;
    auto& sh = L.sheet_parts[0][idx];
    uint32_t rotation_i = sh.get_rotation();
    auto rotate = geo::get_rotate(rotation_i, sh.allowed_rotations);
    Polygon_with_holes_2 rotated = rotate ? geo::transform_polygon_with_holes(*rotate, *sh.base) : *sh.base;
    auto first_v = rotated.outer_boundary().vertices_begin();
    double max_sq_d = 0.0;
    for (auto v = rotated.outer_boundary().vertices_begin(); v != rotated.outer_boundary().vertices_end(); ++v) {
        double dx = CGAL::to_double(v->x() - first_v->x());
        double dy = CGAL::to_double(v->y() - first_v->y());
        double d2 = dx*dx + dy*dy; if (d2 > max_sq_d) max_sq_d = d2;
    }
    double rB = std::sqrt(max_sq_d);
    geo::FT Rft = L.sheets[0].get_radius(); double R = CGAL::to_double(Rft);
    double innerR = R - rB; if (innerR <= 0.0) return;
    const int RINGS = 6; const int ANGULAR = 32;
    double best_tx = 0.0, best_ty = 0.0; double best_pd = std::numeric_limits<double>::infinity(); bool placed = false;
    for (int ring = 0; ring <= RINGS && !placed; ++ring) {
        double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);
        int samples = (ring == 0) ? 1 : ANGULAR;
        for (int s = 0; s < samples; ++s) {
            double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);
            double centre_x = R + ring_r * std::cos(theta);
            double centre_y = R + ring_r * std::sin(theta);
            double tx = centre_x - CGAL::to_double(first_v->x());
            double ty = centre_y - CGAL::to_double(first_v->y());
            sh.set_translate(tx, ty);
            if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;
            double total_pd = 0.0;
            for (size_t j = 0; j < L.sheet_parts[0].size(); ++j) {
                if (j == idx) continue;
                double ox = L.sheet_parts[0][j].get_translate_double_x(); double oy = L.sheet_parts[0][j].get_translate_double_y();
                if (ox < -1e5 || oy < -1e5) continue;
                auto& nfp = comp_nfp(L.sheet_parts[0][j].base, L.sheet_parts[0][j].get_rotation(), L.sheet_parts[0][j].allowed_rotations,
                                      sh.base, rotation_i, sh.allowed_rotations, L);
                FT pd_ft = comp_pd(nfp, tx - ox, ty - oy, L);
                double pd_d = CGAL::to_double(pd_ft);
                total_pd += pd_d; if (total_pd > best_pd) break;
            }
            if (total_pd <= geo::BIAS) { best_tx = tx; best_ty = ty; best_pd = total_pd; placed = true; break; }
            if (total_pd < best_pd) { best_pd = total_pd; best_tx = tx; best_ty = ty; }
        }
    }
    if (placed) {
        sh.set_translate(best_tx, best_ty);
        L.mark_part_placed(idx);
        double placed_x = CGAL::to_double(sh.get_translate_ft_x());
        double placed_y = CGAL::to_double(sh.get_translate_ft_y());
        for (size_t j = 0; j < L.sheet_parts[0].size(); ++j) {
            if (j == idx) continue;
            double ox = L.sheet_parts[0][j].get_translate_double_x(); double oy = L.sheet_parts[0][j].get_translate_double_y();
            auto& nfp = comp_nfp(L.sheet_parts[0][j].base, L.sheet_parts[0][j].get_rotation(), L.sheet_parts[0][j].allowed_rotations,
                                  sh.base, sh.get_rotation(), sh.allowed_rotations, L);
            FT pd = comp_pd(nfp, placed_x - ox, placed_y - oy, L);
            L.set_pd(idx, j, pd);
        }
    } else {
        sh.set_translate(-1e6, -1e6);
        L.mark_part_unplaced(idx);
    }
}

bool minimize_overlap(Layout& layout, volatile bool* requestQuit) {
    if (layout.sheet_parts.empty() || layout.poly_num == 0) return false;
    size_t iter = 0;
    while (iter < layout.maxIterations) {
        if (*requestQuit) return false;
        bool any_improved = false;
        for (size_t i = 0; i < layout.poly_num; ++i) {
            double cur_pd = CGAL::to_double(layout.get_one_polygon_pd(i));
            if (cur_pd <= geo::BIAS) continue;
            auto& sh = layout.sheet_parts[0][i];
            const int TRIES = 20;
            double best_pd = cur_pd;
            double best_tx = sh.get_translate_double_x();
            double best_ty = sh.get_translate_double_y();

            if (layout.sheets[0].is_circle()) {
                geo::FT Rft = layout.sheets[0].get_radius(); double R = CGAL::to_double(Rft);
                auto first_v = sh.base->outer_boundary().vertices_begin();
                double max_sq = 0.0;
                for (auto v = sh.base->outer_boundary().vertices_begin(); v != sh.base->outer_boundary().vertices_end(); ++v) {
                    double dx = CGAL::to_double(v->x() - first_v->x()); double dy = CGAL::to_double(v->y() - first_v->y());
                    double d2 = dx*dx + dy*dy; if (d2 > max_sq) max_sq = d2;
                }
                double rB = std::sqrt(max_sq);
                double innerR = R - rB;
                if (innerR <= 0) {
                    sh.set_translate(-1e6, -1e6);
                    layout.mark_part_unplaced(i);
                    any_improved = true;
                    continue;
                }
                for (int t = 0; t < TRIES; ++t) {
                    double ru = (double)std::rand() / RAND_MAX;
                    double rr = innerR * std::sqrt(ru);
                    double theta = 2.0 * PI * ((double)std::rand() / RAND_MAX);
                    double cx = R + rr * std::cos(theta);
                    double cy = R + rr * std::sin(theta);
                    double tx = cx - CGAL::to_double(first_v->x());
                    double ty = cy - CGAL::to_double(first_v->y());
                    sh.set_translate(tx, ty);
                    if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;
                    double total_pd = 0.0;
                    for (size_t j = 0; j < layout.poly_num; ++j) {
                        if (i == j) continue;
                        auto& other = layout.sheet_parts[0][j];
                        double ox = other.get_translate_double_x(); double oy = other.get_translate_double_y();
                        if (ox < -1e5 || oy < -1e5) continue;
                        auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                            sh.base, sh.get_rotation(), sh.allowed_rotations, layout);
                        FT pd = comp_pd(nfp, tx - ox, ty - oy, layout);
                        total_pd += CGAL::to_double(pd);
                        if (total_pd >= best_pd) break;
                    }
                    if (total_pd < best_pd) { best_pd = total_pd; best_tx = tx; best_ty = ty; }
                }
            } else {
                double height = CGAL::to_double(layout.sheets[0].get_height());
                double length = CGAL::to_double(layout.cur_length);
                if (length <= 0 || height <= 0) continue;
                for (int t = 0; t < TRIES; ++t) {
                    double tx = ((double)std::rand() / RAND_MAX) * length;
                    double ty = ((double)std::rand() / RAND_MAX) * height;
                    sh.set_translate(tx, ty);
                    double total_pd = 0.0;
                    for (size_t j = 0; j < layout.poly_num; ++j) {
                        if (i == j) continue;
                        auto& other = layout.sheet_parts[0][j];
                        double ox = other.get_translate_double_x(); double oy = other.get_translate_double_y();
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

            if (best_pd < cur_pd) {
                sh.set_translate(best_tx, best_ty);
                for (size_t j = 0; j < layout.poly_num; ++j) {
                    if (i == j) continue;
                    auto& other = layout.sheet_parts[0][j];
                    double ox = other.get_translate_double_x(); double oy = other.get_translate_double_y();
                    auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                        sh.base, sh.get_rotation(), sh.allowed_rotations, layout);
                    FT pd = comp_pd(nfp, best_tx - ox, best_ty - oy, layout);
                    layout.set_pd(i, j, pd);
                }
                any_improved = true;
            }
        }
        double cur_total = CGAL::to_double(layout.get_pure_total_pd());
        if (cur_total < geo::BIAS) return true;
        if (!any_improved) break;
        ++iter;
    }
    return false;
}

bool try_place_part_greedy(Layout& layout, const Polygon_with_holes_2& part, double /*R*/) {
    if (layout.sheets.empty() || !layout.sheets[0].is_circle()) return false;
    geo::FT Rft = layout.sheets[0].get_radius(); double R = CGAL::to_double(Rft);
    const int RINGS = 6; const int ANGULAR = 32;
    geo::Polygon_2 outer = part.outer_boundary();
    auto first_v = outer.vertices_begin();
    double max_sq = 0.0;
    for (auto v = outer.vertices_begin(); v != outer.vertices_end(); ++v) {
        double dx = CGAL::to_double(v->x() - first_v->x()); double dy = CGAL::to_double(v->y() - first_v->y()); double d2 = dx*dx + dy*dy; if (d2 > max_sq) max_sq = d2;
    }
    double rB = std::sqrt(max_sq); double innerR = R - rB; if (innerR <= 0.0) return false;
    for (int ring = 0; ring <= RINGS; ++ring) {
        double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);
        int samples = (ring == 0) ? 1 : ANGULAR;
        for (int s = 0; s < samples; ++s) {
            double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);
            double centre_x = R + ring_r * std::cos(theta); double centre_y = R + ring_r * std::sin(theta);
            bool inside = true;
            for (auto v = outer.vertices_begin(); v != outer.vertices_end(); ++v) {
                double dx = CGAL::to_double(v->x()) - centre_x;
                double dy = CGAL::to_double(v->y()) - centre_y;
                if (dx*dx + dy*dy > R*R + 1e-9) { inside = false; break; }
            }
            if (inside) return true;
        }
    }
    return false;
}

static double compute_current_utilization(const Layout& layout) {
    if (layout.sheets.empty() || !layout.sheets[0].is_circle()) return 0.0;
    geo::FT R = layout.sheets[0].get_radius(); double parts_area = 0.0;
    for (size_t i = 0; i < layout.poly_num; ++i) {
        const auto& p = layout.sheet_parts[0][i];
        if (all_vertices_inside_circle(p.transformed, R)) parts_area += CGAL::to_double(geo::pwh_area(p.transformed));
    }
    double Rd = CGAL::to_double(R);
    double sheet_area = PI * Rd * Rd;
    if (sheet_area <= 0.0) return 0.0;
    double util = parts_area / sheet_area;
    if (util > 1.0) util = 1.0;
    return util;
}

static void perturb_random_part(Layout& layout) {
    if (layout.sheet_parts.empty() || layout.sheet_parts[0].empty()) return;
    size_t n = layout.sheet_parts[0].size();
    size_t idx = (size_t)(std::rand() % n);
    auto& sh = layout.sheet_parts[0][idx];
    double jitter = 1.0 + 10.0 * ((double)std::rand() / (double)RAND_MAX);
    double ang = 2.0 * PI * ((double)std::rand() / (double)RAND_MAX);
    double tx = sh.get_translate_double_x() + jitter * std::cos(ang);
    double ty = sh.get_translate_double_y() + jitter * std::sin(ang);
    sh.set_translate(tx, ty);
}

static void try_add_new_part(Layout& layout) { attempt_place_unplaced_part(layout); }
static void try_remove_small_part(Layout& layout) {
    if (layout.sheet_parts.empty() || layout.sheet_parts[0].empty()) return;
    size_t n = layout.sheet_parts[0].size();
    size_t min_idx = 0; double min_a = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i) {
        double a = CGAL::to_double(geo::pwh_area(layout.sheet_parts[0][i].transformed));
        if (a < min_a) { min_a = a; min_idx = i; }
    }
    const double FAR = -1e6;
    layout.sheet_parts[0][min_idx].set_translate(FAR, FAR);
    layout.mark_part_unplaced(min_idx);
}

static void simulated_annealing_optimization(Layout& layout,
    size_t max_time,
    const clock_t& start,
    std::function<void(const Solution&)> ProgressHandler,
    volatile bool* requestQuit) {
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
    auto swap_two_parts = [&](Layout& L) {
        if (L.sheet_parts.empty() || L.sheet_parts[0].size() < 2) return;
        size_t n = L.sheet_parts[0].size();
        size_t a = (size_t)std::min<double>(std::floor(random_double() * n), (double)n - 1);
        size_t b = (size_t)std::min<double>(std::floor(random_double() * n), (double)n - 1);
        if (a == b) return;
        auto& A = L.sheet_parts[0][a];
        auto& B = L.sheet_parts[0][b];
        double ax = A.get_translate_double_x(); double ay = A.get_translate_double_y();
        double bx = B.get_translate_double_x(); double by = B.get_translate_double_y();
        A.set_translate(bx, by);
        B.set_translate(ax, ay);
    };

    while (((double)(clock() - start) / CLOCKS_PER_SEC) < max_time) {
        if (*requestQuit) return;
        for (int i = 0; i < 100; ++i) {
            if (*requestQuit) return;
            Layout new_layout = current_layout;
            int operation = (int)std::floor(random_double() * 5.0);
            switch (operation) {
                case 0: perturb_random_part(new_layout); break;
                case 1: swap_two_parts(new_layout); break;
                case 2: try_add_new_part(new_layout); break;
                case 3: try_remove_small_part(new_layout); break;
                case 4: attempt_place_unplaced_part(new_layout); break;
                default: perturb_random_part(new_layout); break;
            }
            double new_util = compute_current_utilization(new_layout);
            current_util = compute_current_utilization(current_layout);
            double delta_energy = current_util - new_util;
            bool accept = false;
            if (delta_energy < 0.0) accept = true;
            else if (temperature > 0.0) {
                double prob = std::exp(-delta_energy / temperature);
                if (prob > random_double()) accept = true;
            }
            if (accept) {
                current_layout = new_layout;
                current_util = new_util;
                if (new_util > best_util) {
                    best_util = new_util;
                    best_layout = new_layout;
                    double disp_len = 0.0;
                    if (!best_layout.sheets.empty()) disp_len = CGAL::to_double(best_layout.sheets[0].get_diameter());
                    ProgressHandler(Solution(disp_len, best_util, ((double)(clock() - start) / CLOCKS_PER_SEC), best_layout.sheet_parts[0]));
                }
            }
        }
        temperature *= cooling_rate;
        if (temperature < final_temperature) temperature = initial_temperature * 0.5;
        if (((double)(clock() - start) / CLOCKS_PER_SEC) > max_time * 0.8) break;
    }
    layout = best_layout;
    layout.best_utilization = best_util;
}

void GOMH(Layout& layout,
    size_t max_time,
    std::function<void(const Solution&)> ProgressHandler,
    volatile bool* requestQuit) {
    clock_t start = clock();
    double fixed_radius = 150.0;
    if (!layout.sheets.empty()) layout.sheets[0].set_radius(geo::FT(fixed_radius));
    get_greedy_initial_solution(layout);
    double length_display = CGAL::to_double(layout.best_length);
    if (!layout.sheets.empty() && layout.sheets[0].is_circle()) {
        length_display = CGAL::to_double(layout.sheets[0].get_diameter());
    }
    double current_util = compute_current_utilization(layout);
    layout.best_utilization = current_util;
    layout.best_result = layout.sheet_parts;
    ProgressHandler(Solution(length_display, layout.best_utilization, ((double)(clock() - start) / CLOCKS_PER_SEC), layout.best_result[0]));
    shrink(layout);
    while (((double)(clock() - start) / CLOCKS_PER_SEC) < max_time) {
        if (*requestQuit) return;
        auto feasible = minimize_overlap(layout, requestQuit);
        if (feasible) {
            layout.best_result = layout.sheet_parts;
            layout.best_utilization = compute_util(layout, 0.0);
            double disp_len = length_display;
            ProgressHandler(Solution(disp_len, layout.best_utilization, ((double)(clock() - start) / CLOCKS_PER_SEC), layout.best_result[0]));
            shrink(layout);
        } else {
            shrink(layout);
        }
    }
    simulated_annealing_optimization(layout, max_time, start, ProgressHandler, requestQuit);
}

static void get_greedy_initial_solution(Layout& layout) {
    if (layout.sheets.empty() || !layout.sheets[0].is_circle()) throw std::runtime_error("get_greedy_initial_solution requires a circular sheet");
    geo::FT Rft = layout.sheets[0].get_radius(); double R = CGAL::to_double(Rft); double center_x = R, center_y = R;
    size_t total_parts = layout.sheet_parts.size() ? layout.sheet_parts[0].size() : 0;
    if (total_parts == 0) return;
    struct PartArea { size_t idx; double area; };
    std::vector<PartArea> parts; parts.reserve(total_parts);
    for (size_t i = 0; i < total_parts; ++i) parts.push_back({ i, CGAL::to_double(geo::pwh_area(layout.sheet_parts[0][i].transformed)) });
    std::sort(parts.begin(), parts.end(), [](const PartArea& A, const PartArea& B) { return A.area > B.area; });
    std::vector<char> placed(total_parts, 0);
    const int RINGS = 8;
    const int ANGULAR = 64;
    for (const auto& pa : parts) {
        size_t i = pa.idx;
        auto& sh = layout.sheet_parts[0][i];
        uint32_t rotation_i = sh.get_rotation();
        auto rotate = geo::get_rotate(rotation_i, sh.allowed_rotations);
        Polygon_with_holes_2 rotated = rotate ? geo::transform_polygon_with_holes(*rotate, *sh.base) : *sh.base;
        auto first_v = rotated.outer_boundary().vertices_begin();
        double max_sq = 0.0;
        for (auto vit = rotated.outer_boundary().vertices_begin(); vit != rotated.outer_boundary().vertices_end(); ++vit) {
            double dx = CGAL::to_double(vit->x() - first_v->x());
            double dy = CGAL::to_double(vit->y() - first_v->y());
            double d2 = dx * dx + dy * dy;
            if (d2 > max_sq) max_sq = d2;
        }
        double rB = std::sqrt(max_sq);
        double innerR = R - rB;
        if (innerR <= 0.0) continue;
        bool success = false;
        for (int ring = 0; ring <= RINGS && !success; ++ring) {
            double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);
            int samples = (ring == 0) ? 1 : ANGULAR;
            for (int s = 0; s < samples; ++s) {
                double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);
                double cx = center_x + ring_r * std::cos(theta);
                double cy = center_y + ring_r * std::sin(theta);
                double tx = cx - CGAL::to_double(first_v->x());
                double ty = cy - CGAL::to_double(first_v->y());
                sh.set_translate(tx, ty);
                if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;
                double overlap_sum = 0.0;
                for (size_t j = 0; j < total_parts; ++j) {
                    if (!placed[j]) continue;
                    auto& other = layout.sheet_parts[0][j];
                    double ox = other.get_translate_double_x(); double oy = other.get_translate_double_y();
                    auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                         sh.base, rotation_i, sh.allowed_rotations, layout);
                    FT pd = comp_pd(nfp, tx - ox, ty - oy, layout);
                    overlap_sum += CGAL::to_double(pd);
                    if (overlap_sum > geo::BIAS) break;
                }
                if (overlap_sum <= geo::BIAS) { placed[i] = 1; success = true; break; }
            }
        }
    }
    const double FAR = -1e6;
    for (size_t i = 0; i < total_parts; ++i) if (!placed[i]) layout.sheet_parts[0][i].set_translate(FAR, FAR);
    for (size_t i = 0; i < total_parts; ++i) {
        if (!placed[i]) continue;
        for (size_t j = i + 1; j < total_parts; ++j) {
            if (!placed[j]) continue;
            auto& A = layout.sheet_parts[0][i];
            auto& B = layout.sheet_parts[0][j];
            double ax = A.get_translate_double_x(); double ay = A.get_translate_double_y();
            double bx = B.get_translate_double_x(); double by = B.get_translate_double_y();
            auto& nfp = comp_nfp(A.base, A.get_rotation(), A.allowed_rotations, B.base, B.get_rotation(), B.allowed_rotations, layout);
            FT pd = comp_pd(nfp, ax - bx, ay - by, layout);
            layout.set_pd(i, j, pd);
        }
    }
    size_t count = 0; double totalArea = 0.0;
    for (size_t i = 0; i < total_parts; ++i) {
        if (placed[i]) { ++count; totalArea += CGAL::to_double(geo::pwh_area(layout.sheet_parts[0][i].transformed)); }
    }
    std::clog << "Greedy placed " << count << " of " << total_parts << " parts, total area " << totalArea << std::endl;
    layout.update_cur_length();
    layout.best_length = layout.cur_length;
    layout.best_result = layout.sheet_parts;
    layout.best_utilization = compute_util(layout, 0.0);
}

} // namespace nesting
