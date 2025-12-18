#include "nesting.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <thread>
#include <cstdlib>

// Define PI for MSVC where M_PI may be unavailable
#ifndef PI
static constexpr double PI = 3.14159265358979323846;
#endif

namespace nesting {
    // Compute utilization for circular sheets: parts total area / circle area. Cap at 1.0.
    // Keep signature unchanged; for non-circle sheets falls back to rectangular computation.
    static double compute_util(const Layout& layout, double /*length*/) {
        if (layout.sheets.empty()) return 0.0;
        // If circular sheet, compute utilization as parts total area / circle area
        if (layout.sheets[0].is_circle()) {
            double parts_area = CGAL::to_double(layout.area);
            double R = CGAL::to_double(layout.sheets[0].get_radius());
            if (R <= 0.0) return 0.0;
            double sheet_area = PI * R * R;
            if (sheet_area <= 0.0) return 0.0;
            double util = parts_area / sheet_area;
            if (util > 1.0) util = 1.0;
            return util;
        }
        // Fallback for rectangular sheets: use current length * height
        double area_d = CGAL::to_double(layout.area);
        double height = CGAL::to_double(layout.sheets[0].get_height());
        double length = CGAL::to_double(layout.cur_length);
        if (length <= 0.0 || height <= 0.0) return 0.0;
        const double EPS = 1e-12;
        double rect_area = length * height + EPS;
        double util = area_d / rect_area;
        if (util > 1.0) util = 1.0;
        return util;
    }

    // Helper: check all vertices of a polygon with holes lie inside circle centered at (R,R)
    static bool all_vertices_inside_circle(const Polygon_with_holes_2& pwh, const geo::FT& R) {
        if (R <= geo::FT(0)) return false;
        double Rd = CGAL::to_double(R);
        double R2 = Rd * Rd;
        // check outer boundary
        for (auto vit = pwh.outer_boundary().vertices_begin(); vit != pwh.outer_boundary().vertices_end(); ++vit) {
            double dx = CGAL::to_double(vit->x()) - Rd;
            double dy = CGAL::to_double(vit->y()) - Rd;
            if (dx*dx + dy*dy > R2 + 1e-12) return false;
        }
        // check holes
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
        Polygon_with_holes_2 nfp;
        auto kv = layout.nfp_cache.find(key);
        if (kv == layout.nfp_cache.end()) {
            Transformation scale(CGAL::SCALING, -1);
            auto rotate_A = geo::get_rotate(rotation_A, allowed_rotation_A);
            auto rotate_B = geo::get_rotate(rotation_B, allowed_rotation_B);
            Polygon_with_holes_2 minus_B;
            Polygon_with_holes_2 rotated_A;
            if (rotate_B) {
                minus_B = geo::transform_polygon_with_holes(scale * (*rotate_B), *poly_B);
            }
            else {
                minus_B = geo::transform_polygon_with_holes(scale, *poly_B);
            }
            if (rotate_A) {
                rotated_A = geo::transform_polygon_with_holes(*rotate_A, *poly_A);
            }
            else {
                rotated_A = *poly_A;
            }
            nfp = CGAL::minkowski_sum_2(rotated_A, minus_B);
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
        else {
            return kv->second;
        }
    }

    FT comp_pd(NFPCacheValue& v, double px, double py, Layout& layout) {
        auto& original_nfp = v.nfp;
        if (px <= v.xmin || px >= v.xmax || py <= v.ymin || py >= v.ymax) {
            return geo::FT_ZERO;
        }
        hash::PDCacheKey key(&original_nfp, px, py);
        FT pd;
        layout.pd_count++;
        auto iter = layout.pd_cache.find(key);
        if (iter != layout.pd_cache.end()) {
            pd = iter->second;
        }
        else {
            layout.pd_miss++;
            // 为效率考虑使用不精确构建
            Point_2 relative_point(px, py);
            pd = geo::comp_pd(original_nfp, relative_point);
            layout.pd_cache.emplace(key, pd);
        }
        return pd;
    }

    void shrink(Layout& layout) {
        // For circular sheets, re-place any polygons whose vertices fall outside the circle.
        if (layout.sheets.empty()) return;
        if (!layout.sheets[0].is_circle()) {
            // For non-circular (rectangular) sheets we do not perform any length-based shrink here
            return;
        }
        std::vector<size_t> moved_indices;
        // radius and center at (R,R)
        geo::FT R = layout.sheets[0].get_radius();
        double sheetR_d = CGAL::to_double(R);
        for (size_t i = 0; i < layout.poly_num; i++) {
            auto& p = layout.sheet_parts[0][i];
            if (!all_vertices_inside_circle(p.transformed, R)) {
                // compute IFR similarly to previous circular logic
                geo::Polygon_2 ifr_poly;
                auto first_v = p.transformed.outer_boundary().vertices_begin();
                double max_sq_d = 0.0;
                for (auto vtx = p.transformed.outer_boundary().vertices_begin();
                    vtx != p.transformed.outer_boundary().vertices_end(); ++vtx) {
                    double dx = CGAL::to_double(vtx->x() - first_v->x());
                    double dy = CGAL::to_double(vtx->y() - first_v->y());
                    double d2 = dx * dx + dy * dy;
                    if (d2 > max_sq_d) max_sq_d = d2;
                }
                geo::FT rB = geo::FT(std::sqrt(max_sq_d));
                geo::FT innerR = R - rB;
                if (innerR <= geo::FT(0)) {
                    ifr_poly.clear();
                } else {
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
                    Transformation translate(CGAL::TRANSLATION, Vector_2(-first_v->x(), -first_v->y()));
                    ifr_poly = transform(translate, ifr_poly);
                }
                if (ifr_poly.size() > 0) {
                    auto bbox = ifr_poly.bbox();
                    double random_x = bbox.x_span() * rand::right_nd01() + bbox.xmin();
                    double random_y = bbox.y_span() * rand::center_nd01() + bbox.ymin();
                    p.set_translate(random_x, random_y);
                } else {
                    // fallback: place at center
                    p.set_translate(CGAL::to_double(R), CGAL::to_double(R));
                }
                moved_indices.push_back(i);
            }
        }
        // update pd entries for moved polygons
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
        // Update utilization for circular sheet
        layout.best_utilization = compute_util(layout, 0.0);
    }

    // New circular-focused initial placement
    void get_init_solu(Layout& layout) {
        // Circular-sheet initial placement: inner-to-outer concentric-ring sampling
        size_t num_poly = layout.poly_num;
        if (num_poly == 0) return;

        if (layout.sheets.empty() || !layout.sheets[0].is_circle()) {
            throw std::runtime_error("get_init_solu requires a circular sheet");
        }

        // Sort parts by descending area (tie-breaker: max edge length)
        struct Key { size_t idx; double area; double max_edge; };
        std::vector<Key> keys; keys.reserve(num_poly);
        for (size_t i = 0; i < num_poly; ++i) {
            auto& sh = layout.sheet_parts[0][i];
            double a = CGAL::to_double(geo::pwh_area(sh.transformed));
            double max_e = 0.0;
            auto& outer = sh.transformed.outer_boundary();
            for (auto v = outer.vertices_begin(); v != outer.vertices_end(); ++v) {
                auto next = v; ++next;
                if (next == outer.vertices_end()) next = outer.vertices_begin();
                double dx = CGAL::to_double(next->x() - v->x());
                double dy = CGAL::to_double(next->y() - v->y());
                double len = std::sqrt(dx*dx + dy*dy);
                if (len > max_e) max_e = len;
            }
            keys.push_back({i, a, max_e});
        }
        std::sort(keys.begin(), keys.end(), [](const Key& A, const Key& B) {
            if (A.area != B.area) return A.area > B.area; return A.max_edge > B.max_edge;
        });

        // Sheet center (R,R) and radius
        geo::FT Rft = layout.sheets[0].get_radius();
        double R = CGAL::to_double(Rft);

        // Sampling parameters (inner-to-outer)
        const int RINGS = 8;
        const int ANGULAR = 64;

        // Place parts largest-first using concentric ring sampling inside the true circular IFR
        for (size_t ord = 0; ord < num_poly; ++ord) {
            size_t i = keys[ord].idx;
            auto& sh = layout.sheet_parts[0][i];

            uint32_t rotation_i = sh.get_rotation();
            auto rotate = geo::get_rotate(rotation_i, sh.allowed_rotations);
            Polygon_with_holes_2 rotated = rotate ? geo::transform_polygon_with_holes(*rotate, *sh.base) : *sh.base;

            // reference vertex (first) is used as reference point
            auto first_v = rotated.outer_boundary().vertices_begin();

            // compute rB = max distance from reference to polygon vertices
            double max_sq_d = 0.0;
            for (auto v = rotated.outer_boundary().vertices_begin(); v != rotated.outer_boundary().vertices_end(); ++v) {
                double dx = CGAL::to_double(v->x() - first_v->x());
                double dy = CGAL::to_double(v->y() - first_v->y());
                double d2 = dx*dx + dy*dy;
                if (d2 > max_sq_d) max_sq_d = d2;
            }
            double rB = std::sqrt(max_sq_d);

            // innerR: radius available for placing reference point (so whole polygon stays inside sheet)
            double innerR = R - rB;

            bool placed = false;
            double best_tx = 0.0, best_ty = 0.0;
            double best_pd = std::numeric_limits<double>::infinity();

            // Sample from inner rings toward outer rings (center-first)
            for (size_t ring = 0; ring <= (size_t)RINGS && !placed; ++ring) {
                double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);
                int samples = (ring == 0) ? 1 : ANGULAR;
                for (int s = 0; s < samples; ++s) {
                    double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);
                    double centre_x = R + ring_r * std::cos(theta);
                    double centre_y = R + ring_r * std::sin(theta);
                    // translation to move reference vertex to centre_x,centre_y
                    double tx = centre_x - CGAL::to_double(first_v->x());
                    double ty = centre_y - CGAL::to_double(first_v->y());

                    // tentative placement
                    sh.set_translate(tx, ty);

                    // boundary check: every vertex must lie inside the sheet circle
                    if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;

                    // evaluate pd with already placed parts
                    double total_pd = 0.0;
                    for (size_t k = 0; k < ord; ++k) {
                        size_t j = keys[k].idx;
                        auto& other = layout.sheet_parts[0][j];
                        double ox = other.get_translate_double_x();
                        double oy = other.get_translate_double_y();
                        auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                              sh.base, rotation_i, sh.allowed_rotations, layout);
                        FT pd_ft = comp_pd(nfp, tx - ox, ty - oy, layout);
                        double pd_d = CGAL::to_double(pd_ft);
                        total_pd += pd_d;
                        if (total_pd > best_pd) break; // prune bad candidates
                    }

                    // accept perfect candidate immediately
                    if (total_pd <= geo::BIAS) {
                        best_tx = tx; best_ty = ty; best_pd = total_pd; placed = true; break;
                    }
                    // keep best (minimum overlap) candidate
                    if (total_pd < best_pd) { best_pd = total_pd; best_tx = tx; best_ty = ty; }
                }
            }

            // fallback: try placing at sheet center (reference at R,R) then small jitter attempts
            if (!placed) {
                double center_tx = R - CGAL::to_double(first_v->x());
                double center_ty = R - CGAL::to_double(first_v->y());
                sh.set_translate(center_tx, center_ty);
                if (all_vertices_inside_circle(sh.transformed, Rft)) { best_tx = center_tx; best_ty = center_ty; placed = true; }
                else {
                    const int JITTERS = 24;
                    const double JITTER_R = std::max(1.0, std::min(innerR, 5.0));
                    for (int a = 0; a < JITTERS; ++a) {
                        double ang = 2.0 * PI * (double)a / (double)JITTERS;
                        double rx = center_tx + JITTER_R * std::cos(ang);
                        double ry = center_ty + JITTER_R * std::sin(ang);
                        sh.set_translate(rx, ry);
                        if (all_vertices_inside_circle(sh.transformed, Rft)) { best_tx = rx; best_ty = ry; placed = true; break; }
                    }
                    // if still nothing, use best sampled candidate even if overlapping
                    if (!placed && best_pd < std::numeric_limits<double>::infinity()) { sh.set_translate(best_tx, best_ty); placed = true; }
                    // final fallback: leave at center (may overlap)
                    if (!placed) { sh.set_translate(center_tx, center_ty); best_tx = center_tx; best_ty = center_ty; placed = true; }
                }
            } else {
                sh.set_translate(best_tx, best_ty);
            }

            // update PD entries between this placed part and previously placed ones
            double placed_x = CGAL::to_double(sh.get_translate_ft_x());
            double placed_y = CGAL::to_double(sh.get_translate_ft_y());
            for (size_t k = 0; k < ord; ++k) {
                size_t j = keys[k].idx;
                auto& other = layout.sheet_parts[0][j];
                double ox = other.get_translate_double_x();
                double oy = other.get_translate_double_y();
                auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                    sh.base, sh.get_rotation(), sh.allowed_rotations, layout);
                FT pd = comp_pd(nfp, placed_x - ox, placed_y - oy, layout);
                layout.set_pd(i, j, pd);
            }
        }

        // validate initial solution
        auto pure_overlap = layout.get_pure_total_pd();
        if (pure_overlap > geo::BIAS) {
            throw std::runtime_error("Error get_init_solu(): initial solution is not feasible");
        }

        layout.update_cur_length();
        layout.best_length = layout.cur_length;
        layout.best_result = layout.sheet_parts;
        layout.best_utilization = compute_util(layout, 0.0);
    }

    // Placement attempt focused on a circular sheet using greedy concentric rings strategy
    bool try_place_part_greedy(Layout& layout, const Polygon_with_holes_2& part, double /*R*/) {
        if (layout.sheets.empty() || !layout.sheets[0].is_circle()) return false;

        // Sort parts by descending area (tie-breaker: max edge length)
        struct Key { size_t idx; double area; double max_edge; };
        std::vector<Key> keys; keys.reserve(layout.poly_num);
        for (size_t i = 0; i < layout.poly_num; ++i) {
            auto& sh = layout.sheet_parts[0][i];
            double a = CGAL::to_double(geo::pwh_area(sh.transformed));
            double max_e = 0.0;
            auto& outer = sh.transformed.outer_boundary();
            for (auto v = outer.vertices_begin(); v != outer.vertices_end(); ++v) {
                auto next = v; ++next;
                if (next == outer.vertices_end()) next = outer.vertices_begin();
                double dx = CGAL::to_double(next->x() - v->x());
                double dy = CGAL::to_double(next->y() - v->y());
                double len = std::sqrt(dx*dx + dy*dy);
                if (len > max_e) max_e = len;
            }
            keys.push_back({i, a, max_e});
        }
        std::sort(keys.begin(), keys.end(), [](const Key& A, const Key& B) {
            if (A.area != B.area) return A.area > B.area; return A.max_edge > B.max_edge;
        });

        // Sheet center (R,R) and radius
        geo::FT Rft = layout.sheets[0].get_radius();
        double R = CGAL::to_double(Rft);

        // Sampling parameters
        const int RINGS = 8;
        const int ANGULAR = 64;

        // Attempt placement for each part in order of size
        for (size_t ord = 0; ord < keys.size(); ++ord) {
            size_t i = keys[ord].idx;
            auto& sh = layout.sheet_parts[0][i];

            uint32_t rotation_i = sh.get_rotation();
            auto rotate = geo::get_rotate(rotation_i, sh.allowed_rotations);
            Polygon_with_holes_2 rotated = rotate ? geo::transform_polygon_with_holes(*rotate, *sh.base) : *sh.base;

            // reference vertex (first) is used as reference point
            auto first_v = rotated.outer_boundary().vertices_begin();

            // compute rB = max distance from reference to polygon vertices
            double max_sq_d = 0.0;
            for (auto v = rotated.outer_boundary().vertices_begin(); v != rotated.outer_boundary().vertices_end(); ++v) {
                double dx = CGAL::to_double(v->x() - first_v->x());
                double dy = CGAL::to_double(v->y() - first_v->y());
                double d2 = dx*dx + dy*dy;
                if (d2 > max_sq_d) max_sq_d = d2;
            }
            double rB = std::sqrt(max_sq_d);

            // innerR: radius available for placing reference point (so whole polygon stays inside sheet)
            double innerR = R - rB;
            if (innerR <= 0.0) continue; // cannot place this part at all

            bool placed = false;
            double best_tx = 0.0, best_ty = 0.0;
            double best_pd = std::numeric_limits<double>::infinity();

            // Sample from inner rings toward outer rings (center-first)
            for (int ring = 0; ring <= RINGS && !placed; ++ring) {
                double ring_r = (ring == 0) ? 0.0 : (innerR * (double)ring / (double)RINGS);
                int samples = (ring == 0) ? 1 : ANGULAR;
                for (int s = 0; s < samples; ++s) {
                    double theta = (ring == 0) ? 0.0 : (2.0 * PI * s / samples);
                    double centre_x = R + ring_r * std::cos(theta);
                    double centre_y = R + ring_r * std::sin(theta);
                    // translation to move reference vertex to centre_x,centre_y
                    double tx = centre_x - CGAL::to_double(first_v->x());
                    double ty = centre_y - CGAL::to_double(first_v->y());

                    // tentative placement
                    sh.set_translate(tx, ty);

                    // boundary check: every vertex must lie inside the sheet circle
                    if (!all_vertices_inside_circle(sh.transformed, Rft)) continue;

                    // evaluate pd with already placed parts
                    double total_pd = 0.0;
                    for (size_t k = 0; k < ord; ++k) {
                        size_t j = keys[k].idx;
                        auto& other = layout.sheet_parts[0][j];
                        double ox = other.get_translate_double_x();
                        double oy = other.get_translate_double_y();
                        auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                              sh.base, rotation_i, sh.allowed_rotations, layout);
                        FT pd_ft = comp_pd(nfp, tx - ox, ty - oy, layout);
                        double pd_d = CGAL::to_double(pd_ft);
                        total_pd += pd_d;
                        if (total_pd > best_pd) break; // prune bad candidates
                    }

                    // accept perfect candidate immediately
                    if (total_pd <= geo::BIAS) {
                        best_tx = tx; best_ty = ty; best_pd = total_pd; placed = true; break;
                    }
                    // keep best (minimum overlap) candidate
                    if (total_pd < best_pd) { best_pd = total_pd; best_tx = tx; best_ty = ty; }
                }
            }

            // fallback: try placing at sheet center (reference at R,R) then small jitter attempts
            if (!placed) {
                double center_tx = R - CGAL::to_double(first_v->x());
                double center_ty = R - CGAL::to_double(first_v->y());
                sh.set_translate(center_tx, center_ty);
                if (all_vertices_inside_circle(sh.transformed, Rft)) { best_tx = center_tx; best_ty = center_ty; placed = true; }
                else {
                    const int JITTERS = 24;
                    const double JITTER_R = std::max(1.0, std::min(innerR, 5.0));
                    for (int a = 0; a < JITTERS; ++a) {
                        double ang = 2.0 * PI * (double)a / (double)JITTERS;
                        double rx = center_tx + JITTER_R * std::cos(ang);
                        double ry = center_ty + JITTER_R * std::sin(ang);
                        sh.set_translate(rx, ry);
                        if (all_vertices_inside_circle(sh.transformed, Rft)) { best_tx = rx; best_ty = ry; placed = true; break; }
                    }
                    // if still nothing, use best sampled candidate even if overlapping
                    if (!placed && best_pd < std::numeric_limits<double>::infinity()) { sh.set_translate(best_tx, best_ty); placed = true; }
                    // final fallback: leave at center (may overlap)
                    if (!placed) { sh.set_translate(center_tx, center_ty); best_tx = center_tx; best_ty = center_ty; placed = true; }
                }
            } else {
                sh.set_translate(best_tx, best_ty);
            }

            // update PD entries between this placed part and previously placed ones
            double placed_x = CGAL::to_double(sh.get_translate_ft_x());
            double placed_y = CGAL::to_double(sh.get_translate_ft_y());
            for (size_t k = 0; k < ord; ++k) {
                size_t j = keys[k].idx;
                auto& other = layout.sheet_parts[0][j];
                double ox = other.get_translate_double_x();
                double oy = other.get_translate_double_y();
                auto& nfp = comp_nfp(other.base, other.get_rotation(), other.allowed_rotations,
                                    sh.base, sh.get_rotation(), sh.allowed_rotations, layout);
                FT pd = comp_pd(nfp, placed_x - ox, placed_y - oy, layout);
                layout.set_pd(i, j, pd);
            }
        }

        // validate initial solution
        auto pure_overlap = layout.get_pure_total_pd();
        if (pure_overlap > geo::BIAS) {
            throw std::runtime_error("Error get_init_solu(): initial solution is not feasible");
        }

        layout.update_cur_length();
        layout.best_length = layout.cur_length;
        layout.best_result = layout.sheet_parts;
        layout.best_utilization = compute_util(layout, 0.0);
    }

    // compute utilization of currently placed parts (count only parts fully inside circular sheet)
    static double compute_current_utilization(const Layout& layout) {
        if (layout.sheets.empty() || !layout.sheets[0].is_circle()) return 0.0;
        geo::FT R = layout.sheets[0].get_radius();
        double parts_area = 0.0;
        for (size_t i = 0; i < layout.poly_num; ++i) {
            const auto& p = layout.sheet_parts[0][i];
            if (all_vertices_inside_circle(p.transformed, R)) {
                parts_area += CGAL::to_double(geo::pwh_area(p.transformed));
            }
        }
        double Rd = CGAL::to_double(R);
        double sheet_area = PI * Rd * Rd;
        if (sheet_area <= 0.0) return 0.0;
        double util = parts_area / sheet_area;
        if (util > 1.0) util = 1.0;
        return util;
    }

    // perturb a random part's position by a small jitter
    static void perturb_random_part(Layout& L) {
        if (L.sheet_parts.empty() || L.sheet_parts[0].empty()) return;
        size_t n = L.sheet_parts[0].size();
        size_t idx = (size_t)(std::rand() % n);
        auto& sh = L.sheet_parts[0][idx];
        double jitter = 1.0 + 10.0 * ((double)std::rand() / (double)RAND_MAX);
        double ang = 2.0 * PI * ((double)std::rand() / (double)RAND_MAX);
        double tx = sh.get_translate_double_x() + jitter * std::cos(ang);
        double ty = sh.get_translate_double_y() + jitter * std::sin(ang);
        sh.set_translate(tx, ty);
    }

    // attempt to add a previously unplaced part back into the sheet
    static void try_add_new_part(Layout& L) {
        if (L.sheet_parts.empty() || L.sheet_parts[0].empty()) return;
        double R = 0.0;
        if (!L.sheets.empty() && L.sheets[0].is_circle()) R = CGAL::to_double(L.sheets[0].get_radius());
        size_t n = L.sheet_parts[0].size();
        for (size_t i = 0; i < n; ++i) {
            auto& sh = L.sheet_parts[0][i];
            double tx = sh.get_translate_double_x(); double ty = sh.get_translate_double_y();
            if (tx < -1e5 || ty < -1e5) {
                auto first_v = sh.transformed.outer_boundary().vertices_begin();
                if (first_v == sh.transformed.outer_boundary().vertices_end()) break;
                double center_tx = R - CGAL::to_double(first_v->x());
                double center_ty = R - CGAL::to_double(first_v->y());
                double jitter = 1.0 + 4.0 * ((double)std::rand() / (double)RAND_MAX);
                double ang = 2.0 * PI * ((double)std::rand() / (double)RAND_MAX);
                double candx = center_tx + jitter * std::cos(ang);
                double candy = center_ty + jitter * std::sin(ang);
                sh.set_translate(candx, candy);
                if (!L.sheets.empty() && L.sheets[0].is_circle()) {
                    geo::FT Rft = L.sheets[0].get_radius();
                    if (!all_vertices_inside_circle(sh.transformed, Rft)) {
                        sh.set_translate(-1e6, -1e6);
                    }
                }
                break;
            }
        }
    }

    // remove the smallest part by moving it far away
    static void try_remove_small_part(Layout& L) {
        if (L.sheet_parts.empty() || L.sheet_parts[0].empty()) return;
        size_t n = L.sheet_parts[0].size();
        size_t min_idx = 0; double min_a = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < n; ++i) {
            double a = CGAL::to_double(geo::pwh_area(L.sheet_parts[0][i].transformed));
            if (a < min_a) { min_a = a; min_idx = i; }
        }
        const double FAR = -1e6;
        L.sheet_parts[0][min_idx].set_translate(FAR, FAR);
    }

    static void simulated_annealing_optimization(Layout& layout,
        size_t max_time,
        const clock_t& start,
        std::function<void(const Solution&)> ProgressHandler,
        volatile bool* requestQuit) {
        // Parameters
        double initial_temperature = 1.0;
        double final_temperature = 0.01;
        double cooling_rate = 0.95;

        // RNG
        std::mt19937_64 rng((unsigned long)clock() ^ (unsigned long)std::hash<std::thread::id>()(std::this_thread::get_id()));
        auto random_double = [&](double lo = 0.0, double hi = 1.0) {
            std::uniform_real_distribution<double> dist(lo, hi);
            return dist(rng);
        };

        // State
        Layout current_layout = layout;
        Layout best_layout = layout;
        double best_util = compute_current_utilization(layout);
        double current_util = best_util;
        double temperature = initial_temperature;

        // Neighborhood operations (lambdas capture rng indirectly via random_double)
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

        // SA main loop
        while (((double)(clock() - start) / CLOCKS_PER_SEC) < max_time) {
            if (*requestQuit) return;

            // try multiple neighborhood moves per temperature
            for (int i = 0; i < 100; ++i) {
                if (*requestQuit) return;
                Layout new_layout = current_layout;

                int operation = (int)std::floor(random_double() * 4.0);
                switch (operation) {
                    case 0: perturb_random_part(new_layout); break;
                    case 1: swap_two_parts(new_layout); break;
                    case 2: try_add_new_part(new_layout); break;
                    case 3: try_remove_small_part(new_layout); break;
                    default: perturb_random_part(new_layout); break;
                }

                double new_util = compute_current_utilization(new_layout);
                current_util = compute_current_utilization(current_layout);
                double delta_energy = current_util - new_util; // negative => new is better

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

            // cool down
            temperature *= cooling_rate;
            if (temperature < final_temperature) temperature = initial_temperature * 0.5;

            // reserve time for final processing
            if (((double)(clock() - start) / CLOCKS_PER_SEC) > max_time * 0.8) break;
        }

        // commit best
        layout = best_layout;
        layout.best_utilization = best_util;
    }

    // Modified GOMH: fixed circular sheet radius, no shrink() calls, track total_parts_area
    void GOMH(Layout& layout,
        size_t max_time,
        std::function<void(const Solution&)> ProgressHandler,
        volatile bool* requestQuit) {
        clock_t start = clock();

        // 固定板材半径为150，不再调整板材尺寸
        double fixed_radius = 150.0;
        if (!layout.sheets.empty()) {
            layout.sheets[0].set_radius(geo::FT(fixed_radius));
        }

        // 步骤1：生成初始解（尽量放置更多零件）
        get_init_solu(layout);

        // 计算初始利用率
        double current_util = compute_current_utilization(layout);
        layout.best_utilization = current_util;
        layout.best_result = layout.sheet_parts;

        // 添加 total_parts_area 变量，跟踪已放置零件的总面积（完全在板材内）
        double total_parts_area = 0.0;
        if (!layout.sheets.empty() && layout.sheets[0].is_circle()) {
            geo::FT R = layout.sheets[0].get_radius();
            for (size_t i = 0; i < layout.poly_num; ++i) {
                auto& p = layout.sheet_parts[0][i];
                if (all_vertices_inside_circle(p.transformed, R)) {
                    total_parts_area += CGAL::to_double(geo::pwh_area(p.transformed));
                }
            }
        }

        // 报告初始解
        ProgressHandler(Solution(
            2.0 * fixed_radius,  // 显示直径
            current_util,
            ((double)(clock() - start) / CLOCKS_PER_SEC),
            layout.best_result[0]));

        // 步骤2：使用模拟退火优化布局，尝试放置更多零件
        simulated_annealing_optimization(layout, max_time, start, ProgressHandler, requestQuit);

        std::clog << "Best utilization rate: " << layout.best_utilization << std::endl;
    }

}  // namespace nesting
