#include <cstdint>
#include <variant>
#include <array>
#include <vector>
#include <queue>

/*
 * Globals
*/
constexpr int BOARD_SIZE = 20; // board size is 21 x 21

enum class Direction : uint8_t {
    North = 0,
    East  = 1,
    South = 2,
    West  = 3
};

enum class Enclosure : uint8_t {
    None,
    Square,
    Circle
};

enum class Validity : uint8_t {
    Unknown,
    Valid,
    Invalid
};

enum class InBox : uint8_t {
    Unknown,
    No,
    Yes
};

using ArrowMask = uint8_t;

// Bits: 1<<0 = N, 1<<1 = E, 1<<2 = S, 1<<3 = W
constexpr ArrowMask arrow_bit(Direction d) {
    switch (d) {
        case Direction::North: return 1u << 0;
        case Direction::East:  return 1u << 1;
        case Direction::South: return 1u << 2;
        case Direction::West:  return 1u << 3;
    }
    return 0;
}
/*
 * Classes
*/
struct ArrowCell {
    ArrowMask arrows{};   // any combination of N/E/S/W
    uint8_t arrows_distance{1};
};

struct NumberCell {
    uint8_t value{};      // 0–9
    Enclosure enclosure{Enclosure::None};
};

struct EmptyCell {
    // validity will be derived from Cell.valid, see below
};

using CellContent = std::variant<EmptyCell, ArrowCell, NumberCell>;

struct Cell {
    CellContent content{EmptyCell{}};
    Validity validity{Validity::Unknown};
    InBox in_box = InBox::Unknown;         // puzzle variable

    bool is_arrow() const    { return std::holds_alternative<ArrowCell>(content); }
    bool is_number() const   { return std::holds_alternative<NumberCell>(content); }
    bool is_empty() const    { return std::holds_alternative<EmptyCell>(content); }

    NumberCell* as_number()        { return std::get_if<NumberCell>(&content); }
    const NumberCell* as_number() const { return std::get_if<NumberCell>(&content); }

    ArrowCell* as_arrow()        { return std::get_if<ArrowCell>(&content); }
    const ArrowCell* as_arrow() const { return std::get_if<ArrowCell>(&content); }
  
    static Cell make_arrow(std::initializer_list<Direction> dirs) {
        Cell c;
        ArrowMask mask = 0;
        for (auto d : dirs) mask |= arrow_bit(d);
        c.content = ArrowCell{mask};
        c.validity = Validity::Invalid;
        return c;
    }

    static Cell make_number(uint8_t v, Enclosure e = Enclosure::None) {
        Cell c;
        c.content = NumberCell{v, e};
        c.validity = Validity::Valid;
        return c;
    }

    static Cell make_empty() {
        Cell c;
        c.content = EmptyCell{};
        c.validity = Validity::Unknown;
        return c;
    }
};

constexpr int BoardSize = 20;
using Row   = std::array<Cell, BoardSize>;
using Board = std::array<Row, BoardSize>;

inline Cell& at(Board& b, int row, int col) {
    return b[static_cast<std::size_t>(row)]
            [static_cast<std::size_t>(col)];
}

struct Pos {
    int r;
    int c;
};

struct State {
    Board board;
};

/*
* Helpers
*/
inline bool in_bounds(int r, int c) {
    return 0 <= r && r < BoardSize && 0 <= c && c < BoardSize;
}

bool apply_number_constraints(State& s) {
    bool changed = false;

    for (int r = 0; r < BoardSize; ++r) {
        for (int c = 0; c < BoardSize; ++c) {
            const Cell& cell = s.board[r][c];
            if (!cell.is_number()) continue;

            const auto* num = cell.as_number();
            int n = num->value;

            int fixed_yes = 0;
            int unknown   = 0;
            std::vector<Pos> unknown_cells;

            for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                    int rr = r + dr;
                    int cc = c + dc;
                    if (!in_bounds(rr, cc)) continue;

                    const Cell& nb = s.board[rr][cc];
                    switch (nb.in_box) {
                        case InBox::Yes: fixed_yes++; break;
                        case InBox::Unknown:
                            unknown++;
                            unknown_cells.push_back({rr, cc});
                            break;
                        case InBox::No: break;
                    }
                }
            }

            if (fixed_yes > n) return false;
            if (fixed_yes + unknown < n) return false;

            if (fixed_yes == n) {
                for (auto p : unknown_cells) {
                    Cell& nb = s.board[p.r][p.c];
                    if (nb.in_box == InBox::Unknown) {
                        nb.in_box = InBox::No;
                        changed = true;
                    }
                }
            } else if (fixed_yes + unknown == n) {
                for (auto p : unknown_cells) {
                    Cell& nb = s.board[p.r][p.c];
                    if (nb.in_box == InBox::Unknown) {
                        nb.in_box = InBox::Yes;
                        changed = true;
                    }
                }
            }
        }
    }

    return changed;
}

bool check_arrow_integrity(const State& s, int r, int c, ArrowMask mask, uint8_t dist) {
    static const int DR[4] = {-1, 0, 1, 0};
    static const int DC[4] = { 0, 1, 0,-1};

    // 1. Rogue Yes: all 4 directions
    for (int i = 0; i < 4; ++i) {
        for (int d = 1; d < dist; ++d) {
            int rr = r + DR[i] * d;
            int cc = c + DC[i] * d;
            if (!in_bounds(rr, cc)) break;
            if (s.board[rr][cc].in_box == InBox::Yes) {
                return false;
            }
        }
    }

    // 2. Candidate requirement: arrow directions only
    for (int i = 0; i < 4; ++i) {
        if (!((mask >> i) & 1)) continue; // skip if no arrow in this dir

        bool has_candidate = false;
        for (int d = dist; ; ++d) {
            int rr = r + DR[i] * d;
            int cc = c + DC[i] * d;
            if (!in_bounds(rr, cc)) break;

            auto in_box = s.board[rr][cc].in_box;
            if (in_box != InBox::No) {
                has_candidate = true;
                break;
            }
        }
        if (!has_candidate) return false;
    }

    return true;
}

bool apply_arrow_constraints(State& s) {
    bool changed = false;

    for (int r = 0; r < BoardSize; ++r) {
        for (int c = 0; c < BoardSize; ++c) {
            Cell& cell = s.board[r][c];

            if (!cell.is_arrow()) continue;

            if (cell.in_box == InBox::Unknown) {
                cell.in_box = InBox::No;
                changed = true;
            }

            auto* arrow_cell = cell.as_arrow();
            ArrowMask mask = arrow_cell->arrows;
            uint8_t arrows_distance = arrow_cell->arrows_distance;

            bool no_rogue = check_arrow_integrity(s, r, c, mask, arrows_distance);
            if (!no_rogue) return false;

            int rr;
            int cc;

            int fixed_yes  = 0;
            int unknown    = 0;
            std::vector<Pos> unknown_cells;
            int num_arrows = 0;

            for (int i = 0; i < 4; ++i) {
                if (!((mask >> i) & 1)) continue;

                rr = r;
                cc = c;
                
                switch (static_cast<Direction>(i)) {
                    case Direction::North: rr -= arrows_distance; break;
                    case Direction::South: rr += arrows_distance; break;
                    case Direction::East:  cc += arrows_distance; break;
                    case Direction::West:  cc -= arrows_distance; break;
                }

                if (!in_bounds(rr, cc)) continue;
                Cell& nb = s.board[rr][cc];

                ++num_arrows;   
                switch (nb.in_box) {
                    case InBox::Yes: fixed_yes++; break;
                    case InBox::Unknown:
                        unknown++;
                        unknown_cells.push_back({rr, cc});
                        break;
                    case InBox::No: break;
                }
            }

            for (int i = 0; i < 4; ++i) {    
                if ((mask >> i) & 1) continue; // skip arrow directions
                if (fixed_yes == 0) continue; 

                rr = r;
                cc = c;
                
                switch (static_cast<Direction>(i)) {
                    case Direction::North: rr -= arrows_distance; break;
                    case Direction::South: rr += arrows_distance; break;
                    case Direction::East:  cc += arrows_distance; break;
                    case Direction::West:  cc -= arrows_distance; break;
                }

                if (!in_bounds(rr, cc)) continue;
                Cell& nb = s.board[rr][cc];

                if (nb.in_box == InBox::Yes) return false;

                if (nb.in_box == InBox::Unknown) {
                    nb.in_box = InBox::No;
                    changed = true;
                }
            }

            if (fixed_yes + unknown == 0) arrow_cell->arrows_distance += 1;
        }
    }

    return changed;
}

bool check_connectivity_and_holes(const State& s) {
    static const int DR[4] = {-1, 0, 1, 0};
    static const int DC[4] = { 0, 1, 0,-1};

    int start_r = -1;
    int start_c = -1;
    int yes_count = 0;

    int min_r = BoardSize; 
    int max_r = -1;
    int min_c = BoardSize; 
    int max_c = -1;

    for (int r = 0; r < BoardSize; ++r) {
        for (int c = 0; c < BoardSize; ++c) {
            if (s.board[r][c].in_box == InBox::Yes) {
                yes_count++;
                if (start_r == -1) {
                    start_r = r;
                    start_c = c;
                }
                if (r < min_r) min_r = r;
                if (r > max_r) max_r = r;
                if (c < min_c) min_c = c;
                if (c > max_c) max_c = c;
            }
        }
    }

    if (yes_count == 0) {
        return false;
    }

    std::array<std::array<bool, BoardSize>, BoardSize> vis{};
    std::queue<Pos> q;

    vis[start_r][start_c] = true;
    q.push({start_r, start_c});
    int visited_yes = 0;

    while (!q.empty()) {
        auto [r, c] = q.front();
        q.pop();
        visited_yes++;

        for (int i = 0; i < 4; ++i) {
            int rr = r + DR[i];
            int cc = c + DC[i];
            if (!in_bounds(rr, cc)) continue;
            if (vis[rr][cc]) continue;
            if (s.board[rr][cc].in_box != InBox::Yes) continue;

            vis[rr][cc] = true;
            q.push({rr, cc});
        }
    }

    if (visited_yes != yes_count) {
        return false;
    }

    std::array<std::array<bool, BoardSize>, BoardSize> seen_empty{};
    std::queue<Pos> q2;

    auto is_empty = [&](int r, int c) {
        return s.board[r][c].in_box != InBox::Yes; 
    };

    for (int r = min_r; r <= max_r; ++r) {
        for (int c = min_c; c <= max_c; ++c) {
            bool on_border = (r == min_r || r == max_r || c == min_c || c == max_c);
            if (!on_border) continue;
            if (!is_empty(r, c)) continue;
            if (seen_empty[r][c]) continue;

            seen_empty[r][c] = true;
            q2.push({r, c});

            while (!q2.empty()) {
                auto [cr, cc] = q2.front();
                q2.pop();

                for (int i = 0; i < 4; ++i) {
                    int rr = cr + DR[i];
                    int cc2 = cc + DC[i];
                    if (rr < min_r || rr > max_r || cc2 < min_c || cc2 > max_c) continue;
                    if (seen_empty[rr][cc2]) continue;
                    if (!is_empty(rr, cc2)) continue;

                    seen_empty[rr][cc2] = true;
                    q2.push({rr, cc2});
                }
            }
        }
    }

    for (int r = min_r; r <= max_r; ++r) {
        for (int c = min_c; c <= max_c; ++c) {
            if (is_empty(r, c) && !seen_empty[r][c]) {
                // Found a hole
                return false;
            }
        }
    }

    return true;
}


int main() {
    Board board{};

    // Setting up board
    at(board, 0, 0) = Cell::make_arrow({Direction::East});
    at(board, 0, 4) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 0, 7) = Cell::make_arrow({Direction::East});
    at(board, 0, 11) = Cell::make_arrow({Direction::West});
    at(board, 0, 18) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 1, 4) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 1, 6) = Cell::make_arrow({Direction::East});
    at(board, 1, 8) = Cell::make_number(4, Enclosure::Square);
    at(board, 1, 14) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 2, 1) = Cell::make_arrow({Direction::East});
    at(board, 2, 10) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});
    at(board, 2, 12) = Cell::make_number(5, Enclosure::Circle);
    at(board, 2, 13) = Cell::make_number(5);

    at(board, 3, 5) = Cell::make_arrow({Direction::North});
    at(board, 3, 7) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 3, 16) = Cell::make_arrow({Direction::West});
    at(board, 3, 18) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 4, 2) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 4, 4) = Cell::make_arrow({Direction::North});
    at(board, 4, 6) = Cell::make_number(5);
    at(board, 4, 9) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 4, 14) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 4, 19) = Cell::make_arrow({Direction::West});

    at(board, 5, 8) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});
    at(board, 5, 11) = Cell::make_arrow({Direction::North, Direction::South});
    at(board, 5, 13) = Cell::make_arrow({Direction::South});

    at(board, 6, 0) = Cell::make_arrow({Direction::East});
    at(board, 6, 3) = Cell::make_number(5, Enclosure::Square);
    at(board, 6, 6) = Cell::make_number(6);
    at(board, 6, 9) = Cell::make_number(2);
    at(board, 6, 12) = Cell::make_arrow({Direction::North, Direction::West});
    at(board, 6, 16) = Cell::make_arrow({Direction::North});

    at(board, 7, 2) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 7, 5) = Cell::make_arrow({Direction::North, Direction::East, Direction::South, Direction::West});
    at(board, 7, 18) = Cell::make_arrow({Direction::North});

    at(board, 8, 8) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 8, 10) = Cell::make_arrow({Direction::East, Direction::West});
    at(board, 8, 14) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 8, 17) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 9, 1) = Cell::make_arrow({Direction::North});
    at(board, 9, 3) = Cell::make_arrow({Direction::North, Direction::South});
    at(board, 9, 5) = Cell::make_number(4);
    at(board, 9, 7) = Cell::make_number(7);
    at(board, 9, 18) = Cell::make_number(3, Enclosure::Square);

    at(board, 10, 1) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 10, 12) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 10, 14) = Cell::make_number(5);
    at(board, 10, 16) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 10, 18) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});

    at(board, 11, 2) = Cell::make_number(7, Enclosure::Circle);
    at(board, 11, 5) = Cell::make_arrow({Direction::North, Direction::East, Direction::South});
    at(board, 11, 9) = Cell::make_arrow({Direction::North, Direction::East, Direction::South, Direction::West});
    at(board, 11, 11) = Cell::make_number(5, Enclosure::Square);

    at(board, 12, 1) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 12, 14) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 12, 17) = Cell::make_number(6, Enclosure::Square);

    at(board, 13, 3) = Cell::make_arrow({Direction::South});
    at(board, 13, 7) = Cell::make_number(9);
    at(board, 13, 10) = Cell::make_arrow({Direction::East, Direction::West});
    at(board, 13, 13) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});
    at(board, 13, 16) = Cell::make_arrow({Direction::East});
    at(board, 13, 19) = Cell::make_arrow({Direction::South, Direction::West});

    at(board, 14, 6) = Cell::make_number(4);
    at(board, 14, 8) = Cell::make_number(7, Enclosure::Square);
    at(board, 14, 11) = Cell::make_arrow({Direction::East, Direction::South});

    at(board, 15, 0) = Cell::make_arrow({Direction::East});
    at(board, 15, 5) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 15, 10) = Cell::make_number(4, Enclosure::Circle);
    at(board, 15, 13) = Cell::make_number(7);
    at(board, 15, 15) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 15, 17) = Cell::make_number(4, Enclosure::Circle);

    at(board, 16, 1) = Cell::make_arrow({Direction::South});
    at(board, 16, 3) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 16, 12) = Cell::make_number(7);
    at(board, 16, 14) = Cell::make_number(5, Enclosure::Circle);

    at(board, 17, 6) = Cell::make_arrow({Direction::East});
    at(board, 17, 7) = Cell::make_number(5, Enclosure::Square);
    at(board, 17, 9) = Cell::make_arrow({Direction::East, Direction::South, Direction::West});
    at(board, 17, 18) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});

    at(board, 18, 4) = Cell::make_arrow({Direction::East});
    at(board, 18, 11) = Cell::make_number(4);
    at(board, 18, 13) = Cell::make_arrow({Direction::East, Direction::South, Direction::West});
    at(board, 18, 15) = Cell::make_number(4, Enclosure::Circle);

    at(board, 19, 1) = Cell::make_arrow({Direction::East});
    at(board, 19, 8) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 19, 12) = Cell::make_arrow({Direction::South});
    at(board, 19, 15) = Cell::make_arrow({Direction::South});
    at(board, 19, 19) = Cell::make_arrow({Direction::South, Direction::West});
}

// Objective: Iterate through board cells, solving the validity of each using rules. Once the board has no unknown cells, we must assign each valid cell to one of 6 faces, where traversing forward will either move to another face ("fold", depends on orientation & current face) or stay on the same face
//
// If we traverse cells, our current position and orientation on the board must correspond to some position and orientation on the box
//
// Folding rules: For visualization purposes, say face 1 is the "top", faces 2, 3, 4, & 5 are the N, E, S, and W faces respectively, and 6 is the "bottom" -- each box graph edge corresponds to a face and orientation. A fold-move is when we move forward one cell, given a certain orientation in a cardinal direction of the 2D grid, but choose to "fold" over to another face on the box constructed from the cells
// 1 oriented towards 2 -fold-> 2 oriented towards 6, 1 oriented towards 3 -fold-> 3 oriented towards 6 
// 1 oriented towards 4 -fold-> 4 oriented towards 6, 1 oriented towards 5 -fold-> 5 oriented towards 6 
// 1 -> 6
//
// 2 oriented towards 1 -fold-> 1 oriented towards 4, 2 oriented towards 3 -fold-> 3 oriented towards 4 
// 2 oriented towards 5 -fold-> 5 oriented towards 4, 2 oriented towards 6 -fold-> 6 oriented towards 4 
// 2 -> 4
//
// 3 oriented towards 6 -fold-> 6 oriented towards 5, 2 oriented towards 3 -fold-> 3 oriented towards 5 
// 3 oriented towards 1 -fold-> 1 oriented towards 5, 2 oriented towards 5 -fold-> 5 oriented towards 5 
// 3 -> 5
//
// 4 oriented towards 6 -fold-> 6 oriented towards 2, 2 oriented towards 3 -fold-> 3 oriented towards 2 
// 4 oriented towards 1 -fold-> 1 oriented towards 2, 2 oriented towards 5 -fold-> 5 oriented towards 2 
// 4 -> 2
//
// 5 oriented towards 6 -fold-> 6 oriented towards 3, 2 oriented towards 3 -fold-> 3 oriented towards 3 
// 5 oriented towards 1 -fold-> 1 oriented towards 3, 2 oriented towards 5 -fold-> 5 oriented towards 3 
// 5 -> 3
//
// 6 oriented towards 6 -fold-> 6 oriented towards 1, 2 oriented towards 3 -fold-> 3 oriented towards 1 
// 6 oriented towards 1 -fold-> 1 oriented towards 1, 2 oriented towards 5 -fold-> 5 oriented towards 1 
// 6 -> 1
//
// Evidently, if we are on a cell, assigned to a source face and oriented in some cardinal direction (say a north-facing starting cell is at face 1 (source) oriented towards face 2 (destination)), we can optionally pick a destination face to orient towards by turning in a direction on the grid (right, left, or backwards relative to the current cardinal direction orientation on the grid), then a fold-move will take us to the a position on the destination face oriented towards the source face's opposite face
//
// Say we start at a random valid cell, oriented randomly: assume this combination of position and orientation is assigned to face 1 oriented towards face 2. This is the box's "origin", i.e. face 1 oriented towards face 2 AND after a turn that orients towards either face 4 or face 5 (opposite or left) we MUST FOLD if the move forward is valid. If we want to stay on the face, we may not move forward oriented towards face 4 or face 5 past two orthogonal cell spans that both include the origin cell). From the origin, we can move forward with default orientation (2) to increase the box span in the 2-4 dimension, or move forward oriented towards face 3 to increase the box span in the 3-5 dimension. At a certain width and length, we must choose to "fold". 
//
// Once we fold over a span in both positive and negative directions, we have locked in its length alongside the length of the corresponding (parallel) span on the opposite face. The origin's 4 and 5 oriented edges are folded by default. Given three spans, each with an face-opposite face pair and a "ring" of 4 faces, we check if a ring face has both ends of the span "locked in", then "lock in" the span for the opposite face pair
// After folding, the maximum span for the destination face in the orientation we end up in is the same for a "ring" created by four faces consisting of two opposite face pairs (e.g. 1 -> 2-4, 3-5 -> 6). Every ring has four faces and four edges. Let's say by definition, the origin is on a face in both the X and Y rings, its face and its opposite face are not in the Z ring
//
// When making a folded edge, the "thing" you "lock in" is one of the ends of the span passing through a pair of opposite faces ("ring height"), the end of which is on the face you are folding to
//
// Given these rules and a starting position on the board (with only valid cells), search the space of boxes that are valid given the squares and circles
