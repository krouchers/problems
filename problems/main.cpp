#include <cmath>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <algorithm>

#include "gtest/gtest.h"

const float eps = 2.;

template <typename Iterator>
Iterator remove_duplicates(Iterator first, Iterator last)
{
    auto new_last = first;
    for (auto current = first; current != last; ++current)
    {
        if (std::find(first, new_last, *current) == new_last)
        {
            if (new_last != current)
                *new_last = *current;
            ++new_last;
        }
    }
    return new_last;
}

struct point_t
{
    float x = NAN, y = NAN;
    point_t operator-(const point_t &p) const
    {
        auto res_x = x - p.x;
        auto res_y = y - p.y;
        return {res_x, res_y};
    }
    template <typename T>
    bool operator<(const T p) const
    {
        return x < p && y < p;
    }
    bool is_valid() const
    {
        return (x == x && y == y);
    }
    bool operator==(const point_t &rhs) const
    {
        return x == rhs.x && y == rhs.y;
    }
};

point_t abs(point_t point)
{
    return {std::abs(point.x), std::abs(point.y)};
}

enum class point_positoin
{
    RIGHT,
    LEFT,
};

struct line_t
{
    point_t begin, end;
    float a, b, c;
    line_t(const point_t &first, const point_t &second)
        : begin{first}, end{second},
          a{second.y - first.y},
          b{-(second.x - first.x)},
          c{first.y * second.x - first.x * second.y} {}
    point_t intersection_point(const line_t &arg) const
    {
        float x{(b * arg.c - arg.b * c) / (arg.b * a - b * arg.a)};
        float y{(c * arg.a - arg.c * a) / (arg.b * a - b * arg.a)};
        return {x, y};
    }
    point_positoin position_side_of_point(const point_t p) const
    {
        float P = a * p.x + b * p.y + c;
        return P >= 0 ? point_positoin::RIGHT : point_positoin::LEFT;
    }
    bool operator==(const line_t &rhs) const
    {
        return begin == rhs.begin && end == rhs.end;
    }
};

struct poligon_t
{
    std::vector<line_t> edges_{};
    std::vector<point_t> vertices_{};
    poligon_t(std::initializer_list<point_t> &&vertices)
    {
        auto i = vertices.begin();
        for (; i != vertices.end() - 1; ++i)
        {
            vertices_.push_back(*i);
            edges_.push_back({*i, *(i + 1)});
        }
        vertices_.push_back(*i);
        edges_.push_back({*i, *vertices.begin()});
    }

    poligon_t(const std::vector<point_t> &rhs)
    {
        auto i = rhs.begin();
        for (; i != rhs.end() - 1; ++i)
        {
            vertices_.push_back(*i);
            edges_.push_back({*i, *(i + 1)});
        }
        vertices_.push_back(*i);
        edges_.push_back({*i, *rhs.begin()});
    }

    auto begin() const
    {
        return edges_.begin();
    }
    auto end() const
    {
        return edges_.end();
    }
    bool operator==(const poligon_t &rhs) const
    {
        for (auto &vertice : rhs.vertices_)
        {
            if (std::find_if(vertices_.begin(), vertices_.end(), [&](auto &this_vertice)
                             { return (abs(vertice - this_vertice) < eps && abs(vertice - this_vertice) < eps); }) == vertices_.end())
            {
                return false;
            }
        }
        return true;
    }

    float area() const
    {
        float area{0};
        for (int i = 0; i < vertices_.size() - 1; ++i)
        {
            area += vertices_[i].x * vertices_[i + 1].y;
        }
        area += (*(vertices_.end() - 1)).x * vertices_[0].y;

        for (int i = 0; i < vertices_.size() - 1; ++i)
        {
            area -= vertices_[i].y * vertices_[i + 1].x;
        }
        area -= (*(vertices_.end() - 1)).y * vertices_[0].x;
        return std::abs(area) / 2.0;
    }
};

poligon_t clip_poligon(const poligon_t &poligon, const poligon_t &clipping_area)
{
    std::vector<point_t> intermidiate_step{};
    auto result{poligon};
    for (auto &clipping_edge : clipping_area)
    {
        for (auto &poligon_edge : result)
        {
            auto begin_position{clipping_edge.position_side_of_point(poligon_edge.begin)};
            auto end_position{clipping_edge.position_side_of_point(poligon_edge.end)};
            if (begin_position == point_positoin::LEFT && end_position == point_positoin::RIGHT)
            {
                intermidiate_step.push_back(clipping_edge.intersection_point(poligon_edge));
                intermidiate_step.push_back(poligon_edge.end);
            }
            else if (begin_position == point_positoin::RIGHT && end_position == point_positoin::RIGHT)
            {
                intermidiate_step.push_back(poligon_edge.begin);
                intermidiate_step.push_back(poligon_edge.end);
            }
            else if (begin_position == point_positoin::RIGHT && end_position == point_positoin::LEFT)
            {
                intermidiate_step.push_back(poligon_edge.begin);
                intermidiate_step.push_back(clipping_edge.intersection_point(poligon_edge));
            }
        }
        intermidiate_step.erase(remove_duplicates(intermidiate_step.begin(), intermidiate_step.end()), intermidiate_step.end());
        result = intermidiate_step;
        intermidiate_step.clear();
    }
    return result;
}

int main()
{
    poligon_t poligon{{{100.0, 150.0}, {200.0, 250.0}, {300.0, 200.0}}};
    poligon_t clipping_area{{150, 150}, {150, 200}, {200, 200}, {200, 150}};
    std::cout << clip_poligon(poligon, clipping_area).area();
    testing::InitGoogleTest();
    RUN_ALL_TESTS();
    return 0;
}

TEST(Point, Comparing)
{
    point_t p1{1.0, 1.0};
    bool is_less{p1 < 1.1};
    bool is_not_less{p1 < 0.4};
    EXPECT_EQ(is_less, true);
    EXPECT_EQ(is_not_less, false);
}

TEST(Line, Construction)
{
    line_t line1{{0, 0}, {1, 1}};
    EXPECT_EQ(line1.a, 1);
    EXPECT_EQ(line1.b, -1);
    EXPECT_EQ(line1.c, 0);

    line_t line2{{22, 32}, {12, 14}};
    EXPECT_EQ(line2.a, -18);
    EXPECT_EQ(line2.b, 10);
    EXPECT_EQ(line2.c, 76);
}
TEST(Line, ConstructionVerticalLine)
{
    line_t line{{0, 0}, {0, 1}};
    EXPECT_EQ(line.a, 1);
    EXPECT_EQ(line.b, 0);
    EXPECT_EQ(line.c, 0);
}

TEST(Line, ConstructionGorizontalLine)
{
    line_t line{{0, 0}, {1, 0}};
    EXPECT_EQ(line.a, 0);
    EXPECT_EQ(line.b, -1);
    EXPECT_EQ(line.c, 0);
}

TEST(Poligon, PoligonConstruction)
{
    poligon_t poligon{{{100.0, 150.0}, {200.0, 250.0}, {300.0, 200.0}}};
    EXPECT_EQ(poligon.edges_[0].a, 100);
    EXPECT_EQ(poligon.edges_[0].b, -100);
    EXPECT_EQ(poligon.edges_[0].c, 5000);

    EXPECT_EQ(poligon.edges_[1].a, -50);
    EXPECT_EQ(poligon.edges_[1].b, -100);
    EXPECT_EQ(poligon.edges_[1].c, 35000);

    EXPECT_EQ(poligon.edges_[2].a, -50);
    EXPECT_EQ(poligon.edges_[2].b, 200);
    EXPECT_EQ(poligon.edges_[2].c, -25000);
}

TEST(IntersectionTest, FirstSet)
{
    poligon_t poligon{{{100.0, 150.0}, {200.0, 250.0}, {300.0, 200.0}}};
    poligon_t clipping_area{{150, 150}, {150, 200}, {200, 200}, {200, 150}};
    poligon_t resulting_poligon{{150, 162}, {150, 200}, {200, 200}, {200, 174}};
    auto clipped_poligon = clip_poligon(poligon, clipping_area);
    EXPECT_EQ(resulting_poligon, clipped_poligon);
}

TEST(Line, IntersectionPoint)
{
    poligon_t poligon{{100.0, 150.0}, {200.0, 250.0}, {300.0, 200.0}};
    poligon_t clipping_area{{150, 150}, {150, 200}, {200, 200}, {200, 150}};
    auto intersecton_point = poligon.edges_[0].intersection_point(clipping_area.edges_[0]);
    std::cout << poligon.edges_[0].a << " " << poligon.edges_[0].b << " " << poligon.edges_[0].c << "\n";
    std::cout << clipping_area.edges_[0].a << " " << clipping_area.edges_[0].b << " " << clipping_area.edges_[0].c << "\n";
    point_t expected_value{150, 200};
    EXPECT_EQ(intersecton_point.x, expected_value.x);
    EXPECT_EQ(intersecton_point.y, expected_value.y);
}

TEST(Line, NoIntersectoinPoint)
{
    line_t line1{{0, 0}, {0, 1}};
    line_t line2{{1, 0}, {1, 1}};
    auto ip = line1.intersection_point(line2);
    EXPECT_EQ(ip.is_valid(), false);
}

TEST(Line, PointPosition)
{
    poligon_t poligon{{{100.0, 150.0}, {200.0, 250.0}, {300.0, 200.0}}};
    poligon_t clipping_area{{150, 150}, {150, 200}, {200, 200}, {200, 150}};

    auto positoin1{clipping_area.edges_[0].position_side_of_point(poligon.edges_[0].begin)};
    EXPECT_EQ(positoin1, point_positoin::LEFT);

    auto positoin2{clipping_area.edges_[0].position_side_of_point(poligon.edges_[0].end)};
    EXPECT_EQ(positoin2, point_positoin::RIGHT);
}

TEST(Poligon, Equivalence)
{
    poligon_t poligon1{{100.0, 150.0}, {200.0, 250.0}, {300.0, 200.0}};
    poligon_t poligon2{{200.0, 250.0}, {300.0, 200.0}, {100.0, 150.0}};

    bool testing_value{poligon1 == poligon2};
    EXPECT_EQ(testing_value, true);
}

TEST(Poligon, area)
{
    poligon_t poligon{{-4, 3}, {2, 5}, {5, 1}};
    EXPECT_EQ(poligon.area(), 15);
}
