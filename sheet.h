#pragma once
#include "basic.h"
#include <cmath>
namespace nesting {
    struct Sheet {
        enum class Type { Rectangle, Circle };
        geo::FT width;
        geo::FT height;
        // circle specific
        geo::FT radius;
        Type type{ Type::Circle };
        geo::Polygon_with_holes_2 sheet;
        // 构造函数 - rectangle
        explicit Sheet(geo::FT _width, geo::FT _height)
            : width(_width), height(_height), radius(CGAL::min(_width, _height) / geo::FT(2)), type(Type::Rectangle) {
            build_rectangle_polygon(width, height);
        }
        // 构造函数 - circle
        explicit Sheet(geo::FT _radius)
            : width(geo::FT(2) * _radius), height(geo::FT(2) * _radius), radius(_radius), type(Type::Circle) {
            build_circle_polygon(radius);
        }

        // 将sheet的宽度变为_width (更新为圆时以直径更新为圆)
        inline void set_width(geo::FT _width) {
            if (type == Type::Circle) {
                // treat _width as diameter for circle
                width = _width;
                height = _width;
                radius = _width / geo::FT(2);
                build_circle_polygon(radius);
            }
            else {
                width = _width;
                // rebuild rectangle polygon to reflect new dimensions
                build_rectangle_polygon(width, height);
            }
        }
        // 将sheet的高度变为_height (更新为圆时以直径更新为圆)
        inline void set_height(geo::FT _height) {
            if (type == Type::Circle) {
                // treat _height as diameter for circle
                width = _height;
                height = _height;
                radius = _height / geo::FT(2);
                build_circle_polygon(radius);
            }
            else {
                height = _height;
                build_rectangle_polygon(width, height);
            }
        }

        inline bool is_circle() const { return type == Type::Circle; }

        inline geo::FT area() const {
            if (type == Type::Rectangle) return width * height;
            // circle
            return geo::FT(nesting::geo::PI) * radius * radius;
        }

        inline void set_radius(geo::FT _radius) {
            type = Type::Circle;
            radius = _radius;
            width = geo::FT(2) * radius;
            height = geo::FT(2) * radius;
            build_circle_polygon(radius);
        }

        // Accessors to avoid direct member access elsewhere in the code
        inline geo::FT get_width() const { return width; }
        inline geo::FT get_height() const { return height; }
        inline geo::FT get_radius() const { return radius; }
        inline geo::FT get_diameter() const { return geo::FT(2) * radius; }

    private:
        void build_circle_polygon(geo::FT r) {
            auto& outer = sheet.outer_boundary();
            outer.clear();
            // approximate circle by many-sided polygon
            const int SEGMENTS = 128;
            double rd = CGAL::to_double(r);
            double cx = rd; // center x
            double cy = rd; // center y
            for (int i = 0; i < SEGMENTS; ++i) {
                double theta = 2.0 * geo::PI * static_cast<double>(i) / static_cast<double>(SEGMENTS);
                double x = cx + rd * std::cos(theta);
                double y = cy + rd * std::sin(theta);
                outer.push_back(geo::Point_2(x, y));
            }
            // clear holes
            sheet.holes().clear();
        }

        void build_rectangle_polygon(geo::FT w, geo::FT h) {
            auto& outer = sheet.outer_boundary();
            outer.clear();
            double wd = CGAL::to_double(w);
            double hd = CGAL::to_double(h);
            outer.push_back(geo::Point_2(0, 0));
            outer.push_back(geo::Point_2(wd, 0));
            outer.push_back(geo::Point_2(wd, hd));
            outer.push_back(geo::Point_2(0, hd));
            // clear holes
            sheet.holes().clear();
        }

    };
}  // namespace nesting