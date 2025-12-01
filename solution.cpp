#include <cstdint>
#include <variant>
#include <array>
#include <vector>
#include <queue>
#include <variant>
#include <initializer_list>
#include <iostream>

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

enum class ConstraintStatus { Ok, Contradiction };

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


enum class FaceId : uint8_t {
    F0 = 0, F1, F2, F3, F4, F5
};

struct FaceAssignment {
    bool assigned = false;
    FaceId face{};
    int u = 0;  // local coords on that face
    int v = 0;
};

struct FaceInfo {
    bool used = false;
    int min_u = 0, max_u = -1;
    int min_v = 0, max_v = -1;
    int cell_count = 0;
};

inline bool face_used(uint8_t mask, FaceId f) {
    return mask & (1u << face_index(f));
}

inline uint8_t mark_face_used(uint8_t mask, FaceId f) {
    return mask | (1u << face_index(f));
}

struct FoldingState {
    // For each board cell (r,c), what face / local coord is it on (if any)?
    std::array<std::array<FaceAssignment, BoardSize>, BoardSize> assign{};

    // Info about each face's occupied rectangle
    std::array<FaceInfo, 6> faces{};
};

struct FoldSearchState {
    FoldingState F;
    std::queue<Pos> frontier;
    int yes_total = 0;
    uint8_t used_faces_mask = 0; // bit i set => FaceId::Fi used
};

inline int face_index(FaceId f) {
    return static_cast<int>(f);
}


bool place_cell(FoldingState& F, int r, int c, FaceId face, int u, int v) {
    FaceAssignment& cur = F.assign[r][c];
    int fi = face_index(face);

    // If already assigned, must be consistent
    if (cur.assigned) {
        if (cur.face != face || cur.u != u || cur.v != v) {
            return false; // conflicting assignment
        }
        return true; // nothing new
    }

    // Check overlap: no other cell can already be at (face,u,v)
    for (int rr = 0; rr < BoardSize; ++rr) {
        for (int cc = 0; cc < BoardSize; ++cc) {
            const FaceAssignment& fa = F.assign[rr][cc];
            if (!fa.assigned) continue;
            if (fa.face == face && fa.u == u && fa.v == v) {
                // Another cell already occupies this local coord
                return false;
            }
        }
    }

    // Assign
    cur.assigned = true;
    cur.face = face;
    cur.u = u;
    cur.v = v;

    // Update face info
    FaceInfo& info = F.faces[fi];
    if (!info.used) {
        info.used = true;
        info.min_u = info.max_u = u;
        info.min_v = info.max_v = v;
        info.cell_count = 1;
    } else {
        if (u < info.min_u) info.min_u = u;
        if (u > info.max_u) info.max_u = u;
        if (v < info.min_v) info.min_v = v;
        if (v > info.max_v) info.max_v = v;
        info.cell_count += 1;
    }

    return true;
}

inline bool in_bounds(int r, int c) {
    return 0 <= r && r < BoardSize && 0 <= c && c < BoardSize;
}

bool embed_as_single_face(const State& s, FoldingState& F) {
    F = FoldingState{};

    static const int DR[4] = {-1, 0, 1, 0};
    static const int DC[4] = { 0, 1, 0,-1};

    int start_r = -1, start_c = -1;
    int yes_count = 0;

    for (int r = 0; r < BoardSize; ++r) {
        for (int c = 0; c < BoardSize; ++c) {
            if (s.board[r][c].in_box == InBox::Yes) {
                if (start_r == -1) {
                    start_r = r;
                    start_c = c;
                }
                yes_count++;
            }
        }
    }

    if (yes_count == 0) return false;

    std::array<std::array<bool, BoardSize>, BoardSize> seen{};
    std::queue<Pos> q;

    if (!place_cell(F, start_r, start_c, FaceId::F0, 0, 0)) {
        return false;
    }

    seen[start_r][start_c] = true;
    q.push({start_r, start_c});

    int visited_yes = 0;

    while (!q.empty()) {
        Pos p = q.front();
        q.pop();
        visited_yes++;

        const FaceAssignment& fa = F.assign[p.r][p.c];
        FaceId face = fa.face;
        int u = fa.u;
        int v = fa.v;

        for (int i = 0; i < 4; ++i) {
            int rr = p.r + DR[i];
            int cc = p.c + DC[i];
            if (!in_bounds(rr, cc)) continue;
            if (s.board[rr][cc].in_box != InBox::Yes) continue;
            if (seen[rr][cc]) continue;

            int uu = u + DC[i];
            int vv = v + DR[i];

            if (!place_cell(F, rr, cc, face, uu, vv)) {
                return false;
            }

            seen[rr][cc] = true;
            q.push({rr, cc});
        }
    }

    if (visited_yes != yes_count) {
        return false;
    }

    FaceInfo& info = F.faces[face_index(FaceId::F0)];
    if (!info.used) return false;

    int width  = info.max_u - info.min_u + 1;
    int height = info.max_v - info.min_v + 1;

    if (width * height != info.cell_count) {
        return false;
    }

    return true;
}


bool check_neighbor_consistency(const FoldingState& F,
                                const Pos& p, const Pos& q,
                                int dr, int dc) {
    const auto& ap = F.assign[p.r][p.c];
    const auto& aq = F.assign[q.r][q.c];

    if (!ap.assigned || !aq.assigned) return true; // nothing to check

    if (ap.face == aq.face) {
        int du = aq.u - ap.u;
        int dv = aq.v - ap.v;

        // board move (dr,dc): N=(-1,0),S=(1,0),E=(0,1),W=(0,-1)
        // we map: u = x-like = column, v = y-like = row
        // so dc → du, dr → dv
        if (du == dc && dv == dr && std::abs(du) + std::abs(dv) == 1) {
            return true; // consistent same-face adjacency
        }
        return false; // conflict
    } else {
        // Different faces: for now, don't enforce more.
        return true;
    }
}


static const int DR[4] = {-1, 0, 1, 0};
static const int DC[4] = { 0, 1, 0,-1};

bool final_face_shape_check(const FoldSearchState& S) {
    // For now: just ensure each used face is a solid rectangle
    for (int fi = 0; fi < 6; ++fi) {
        const FaceInfo& info = S.F.faces[fi];
        if (!info.used) continue;
        int w = info.max_u - info.min_u + 1;
        int h = info.max_v - info.min_v + 1;
        if (w * h != info.cell_count) return false;
    }
    // Later: add (A,B,C) side-length compatibility and cube-graph checks
    return true;
}

bool all_yes_assigned(const State& s, const FoldingState& F, int yes_total) {
    int count = 0;
    for (int r = 0; r < BoardSize; ++r)
        for (int c = 0; c < BoardSize; ++c)
            if (s.board[r][c].in_box == InBox::Yes && F.assign[r][c].assigned)
                ++count;
    return count == yes_total;
}

bool search_fold(const State& s, FoldSearchState state) {
    // If no more frontier cells to expand, we are done embedding adjacency
    if (state.frontier.empty()) {
        // All reachable YES cells are now assigned.
        if (!all_yes_assigned(s, state.F, state.yes_total)) {
            return false;
        }
        return final_face_shape_check(state);
    }

    // Pop one cell to process its neighbors
    Pos p = state.frontier.front();
    state.frontier.pop();

    const auto& ap = state.F.assign[p.r][p.c];
    FaceId face_p = ap.face;
    int u_p = ap.u;
    int v_p = ap.v;

    // Explore all neighbors of p
    for (int dir = 0; dir < 4; ++dir) {
        int rr = p.r + DR[dir];
        int cc = p.c + DC[dir];
        if (!in_bounds(rr, cc)) continue;
        if (s.board[rr][cc].in_box != InBox::Yes) continue;

        Pos q{rr, cc};

        const auto& aq = state.F.assign[rr][cc];

        // 1. If neighbor already assigned, just check consistency
        if (aq.assigned) {
            if (!check_neighbor_consistency(state.F, p, q, DR[dir], DC[dir])) {
                return false;
            }
            continue;
        }

        // 2. Neighbor not assigned yet → we branch

        // --- Branch A: put q on the SAME face as p ---
        {
            FoldSearchState branch = state; // copy
            int u_q = u_p + DC[dir];
            int v_q = v_p + DR[dir];

            if (place_cell(branch.F, rr, cc, face_p, u_q, v_q)) {
                branch.frontier.push(q);
                if (search_fold(s, std::move(branch))) {
                    return true;
                }
            }
        }

        // --- Branch B: put q on a NEW face (if any free) ---
        for (int fi = 0; fi < 6; ++fi) {
            FaceId f = static_cast<FaceId>(fi);
            if (face_used(state.used_faces_mask, f)) continue;
            // Don't "move" it to the same face, we already tried that
            if (f == face_p) continue;

            FoldSearchState branch = state; // copy
            branch.used_faces_mask = mark_face_used(branch.used_faces_mask, f);

            // On new face, start q at (0,0) in local coords
            if (place_cell(branch.F, rr, cc, f, 0, 0)) {
                branch.frontier.push(q);
                if (search_fold(s, std::move(branch))) {
                    return true;
                }
            }
        }

        // If no branch returned true for this neighbor, we continue the loop.
        // Important: we *don't* immediately return false; another neighbor
        // or another path might still succeed.
    }

    // After processing all neighbors of p, continue with remaining frontier
    return search_fold(s, std::move(state));
}


bool fold_polyomino_into_faces(const State& s, FoldingState& out_solution) {
    int yes_total = 0;
    Pos start{-1, -1};

    for (int r = 0; r < BoardSize; ++r) {
        for (int c = 0; c < BoardSize; ++c) {
            if (s.board[r][c].in_box == InBox::Yes) {
                yes_total++;
                if (start.r == -1) {
                    start = {r, c};
                }
            }
        }
    }

    if (yes_total == 0) return false;

    FoldSearchState init;
    init.yes_total = yes_total;
    init.used_faces_mask = 0;

    // Place the starting cell on face F0 at (0,0)
    if (!place_cell(init.F, start.r, start.c, FaceId::F0, 0, 0)) {
        return false;
    }
    init.used_faces_mask = mark_face_used(init.used_faces_mask, FaceId::F0);
    init.frontier.push(start);

    if (search_fold(s, init)) {
        out_solution = std::move(init.F);
        return true;
    }

    return false;
}

/*
* Helpers
*/
std::string cell_to_string(const Cell& cell) {
    // 1) Arrow cells: ^ > v < for NESW, possibly in combination
    if (cell.is_arrow()) {
        const ArrowCell* ac = cell.as_arrow();
        ArrowMask mask = ac->arrows;
        std::string s;

        auto has_dir = [&](Direction d) {
            return (mask >> static_cast<int>(d)) & 1u;
        };

        if (has_dir(Direction::North)) s.push_back('^');
        if (has_dir(Direction::East))  s.push_back('>');
        if (has_dir(Direction::South)) s.push_back('v');
        if (has_dir(Direction::West))  s.push_back('<');

        if (s.empty()) s = " "; // shouldn't happen, but just in case
        return s;
    }

    // 2) Number cells: (4), [4], or 4
    if (cell.is_number()) {
        const NumberCell* nc = cell.as_number();
        char d = static_cast<char>('0' + nc->value);
        switch (nc->enclosure) {
            case Enclosure::Circle:
                return std::string{"("} + d + ")";
            case Enclosure::Square:
                return std::string{"["} + d + "]";
            case Enclosure::None:
            default:
                return std::string(1, d);
        }
    }

    // 3) Empty cells: just blank
    return " ";
}

std::string pad_center(const std::string& s, int width) {
    if (static_cast<int>(s.size()) >= width) {
        return s.substr(0, width);
    }

    int total = width - static_cast<int>(s.size());
    int left  = total / 2;
    int right = total - left;

    return std::string(left, ' ') + s + std::string(right, ' ');
}

void print_board(const State& s, std::ostream& os = std::cout) {
    constexpr int CELL_WIDTH = 5; // interior width between '|' and next '|'

    for (int r = 0; r < BoardSize; ++r) {
        for (int c = 0; c < BoardSize; ++c) {
            const Cell& cell = s.board[r][c];
            std::string txt = cell_to_string(cell);

            std::string centered = pad_center(txt, CELL_WIDTH);
            os << '|' << centered;
        }
        os << "|\n";
	for (int c = 0; c < BoardSize; ++c) {
	    os << "------";
	}
	os << "-\n";
    }
}

ConstraintStatus apply_number_constraints(State& s, bool& changed) {
    changed = false;

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

            if (fixed_yes > n) return ConstraintStatus::Contradiction;
            if (fixed_yes + unknown < n) return ConstraintStatus::Contradiction;

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

    return ConstraintStatus::Ok;
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

ConstraintStatus apply_arrow_constraints(State& s, bool& changed) {
    changed = false;

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
            if (!no_rogue) return ConstraintStatus::Contradiction;

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

                if (nb.in_box == InBox::Yes) return ConstraintStatus::Contradiction;

                if (nb.in_box == InBox::Unknown) {
                    nb.in_box = InBox::No;
                    changed = true;
                }
            }

            if (fixed_yes + unknown == 0) arrow_cell->arrows_distance += 1;
        }
    }

    return ConstraintStatus::Ok;
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

    // Get bounding rect
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

bool propagate(State& s) {
    bool any_changed = false;

    while (true) {
        bool changed_this_round = false;

        bool changed_num = false;
        if (apply_number_constraints(s, changed_num) == ConstraintStatus::Contradiction)
            return false;

        bool changed_arrow = false;
        if (apply_arrow_constraints(s, changed_arrow) == ConstraintStatus::Contradiction)
            return false;

        changed_this_round = changed_num || changed_arrow;

        if (!changed_this_round) break;
        any_changed = true;
    }

    return true; // consistent after propagation
}

bool solve(State& s) {
    // 1. Propagate constraints
    if (!propagate(s)) return false;  // contradiction

    // 2. If fully decided
    //!has_unknown_cells(s)
    if (false) {
        // if (!check_connectivity_and_holes(s)) return false;
        // if (!final_exact_arrow_check(s)) return false;
        // later: foldability check
        // success!
        return true;
    }

    // 3. Choose an Unknown cell to branch on
    // Pos p = pick_unknown_cell(s);  // heuristic of your choice

    // Try Yes
    //{
    //    State s_yes = s;
    //    s_yes.board[p.r][p.c].in_box = InBox::Yes;
    //    if (solve(s_yes)) {
    //        s = std::move(s_yes); // keep solution if you want
    //        return true;
    //    }
    //}

    // Try No
    //{
    //    State s_no = s;
    //    s_no.board[p.r][p.c].in_box = InBox::No;
    //    if (solve(s_no)) {
    //        s = std::move(s_no);
    //        return true;
    //    }
    //}

    // Neither branch worked
    return false;
}

int main() {
    Board board{};

    // Setting up board
    at(board, 19, 0) = Cell::make_arrow({Direction::East});
    at(board, 19, 4) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 19, 7) = Cell::make_arrow({Direction::East});
    at(board, 19, 11) = Cell::make_arrow({Direction::West});
    at(board, 19, 18) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 18, 4) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 18, 6) = Cell::make_arrow({Direction::East});
    at(board, 18, 8) = Cell::make_number(4, Enclosure::Square);
    at(board, 18, 14) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 17, 1) = Cell::make_arrow({Direction::East});
    at(board, 17, 10) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});
    at(board, 17, 12) = Cell::make_number(5, Enclosure::Circle);
    at(board, 17, 13) = Cell::make_number(5);

    at(board, 16, 5) = Cell::make_arrow({Direction::North});
    at(board, 16, 7) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 16, 16) = Cell::make_arrow({Direction::West});
    at(board, 16, 18) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 15, 2) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 15, 4) = Cell::make_arrow({Direction::North});
    at(board, 15, 6) = Cell::make_number(5);
    at(board, 15, 9) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 15, 14) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 15, 19) = Cell::make_arrow({Direction::West});

    at(board, 14, 8) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});
    at(board, 14, 11) = Cell::make_arrow({Direction::North, Direction::South});
    at(board, 14, 13) = Cell::make_arrow({Direction::South});

    at(board, 13, 0) = Cell::make_arrow({Direction::East});
    at(board, 13, 3) = Cell::make_number(5, Enclosure::Square);
    at(board, 13, 6) = Cell::make_number(6);
    at(board, 13, 9) = Cell::make_number(2);
    at(board, 13, 12) = Cell::make_arrow({Direction::North, Direction::West});
    at(board, 13, 16) = Cell::make_arrow({Direction::North});

    at(board, 12, 2) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 12, 5) = Cell::make_arrow({Direction::North, Direction::East, Direction::South, Direction::West});
    at(board, 12, 18) = Cell::make_arrow({Direction::North});

    at(board, 11, 8) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 11, 10) = Cell::make_arrow({Direction::East, Direction::West});
    at(board, 11, 14) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 11, 17) = Cell::make_arrow({Direction::North, Direction::West});

    at(board, 10, 1) = Cell::make_arrow({Direction::North});
    at(board, 10, 3) = Cell::make_arrow({Direction::North, Direction::South});
    at(board, 10, 5) = Cell::make_number(4);
    at(board, 10, 7) = Cell::make_number(7);
    at(board, 10, 18) = Cell::make_number(3, Enclosure::Square);

    at(board, 9, 1) = Cell::make_arrow({Direction::North, Direction::East});
    at(board, 9, 12) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 9, 14) = Cell::make_number(5);
    at(board, 9, 16) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 9, 18) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});

    at(board, 8, 2) = Cell::make_number(7, Enclosure::Circle);
    at(board, 8, 5) = Cell::make_arrow({Direction::North, Direction::East, Direction::South});
    at(board, 8, 9) = Cell::make_arrow({Direction::North, Direction::East, Direction::South, Direction::West});
    at(board, 8, 11) = Cell::make_number(5, Enclosure::Square);

    at(board, 7, 1) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 7, 14) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 7, 17) = Cell::make_number(6, Enclosure::Square);

    at(board, 6, 3) = Cell::make_arrow({Direction::South});
    at(board, 6, 7) = Cell::make_number(9);
    at(board, 6, 10) = Cell::make_arrow({Direction::East, Direction::West});
    at(board, 6, 13) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});
    at(board, 6, 16) = Cell::make_arrow({Direction::East});
    at(board, 6, 19) = Cell::make_arrow({Direction::South, Direction::West});

    at(board, 5, 6) = Cell::make_number(4);
    at(board, 5, 8) = Cell::make_number(7, Enclosure::Square);
    at(board, 5, 11) = Cell::make_arrow({Direction::East, Direction::South});

    at(board, 4, 0) = Cell::make_arrow({Direction::East});
    at(board, 4, 5) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 4, 10) = Cell::make_number(4, Enclosure::Circle);
    at(board, 4, 13) = Cell::make_number(7);
    at(board, 4, 15) = Cell::make_arrow({Direction::North, Direction::East, Direction::West});
    at(board, 4, 17) = Cell::make_number(4, Enclosure::Circle);

    at(board, 3, 1) = Cell::make_arrow({Direction::South});
    at(board, 3, 3) = Cell::make_arrow({Direction::East, Direction::South});
    at(board, 3, 12) = Cell::make_number(7);
    at(board, 3, 14) = Cell::make_number(5, Enclosure::Circle);

    at(board, 2, 6) = Cell::make_arrow({Direction::East});
    at(board, 2, 7) = Cell::make_number(5, Enclosure::Square);
    at(board, 2, 9) = Cell::make_arrow({Direction::East, Direction::South, Direction::West});
    at(board, 2, 18) = Cell::make_arrow({Direction::North, Direction::South, Direction::West});

    at(board, 1, 4) = Cell::make_arrow({Direction::East});
    at(board, 1, 11) = Cell::make_number(4);
    at(board, 1, 13) = Cell::make_arrow({Direction::East, Direction::South, Direction::West});
    at(board, 1, 15) = Cell::make_number(4, Enclosure::Circle);

    at(board, 0, 1) = Cell::make_arrow({Direction::East});
    at(board, 0, 8) = Cell::make_arrow({Direction::South, Direction::West});
    at(board, 0, 12) = Cell::make_arrow({Direction::South});
    at(board, 0, 15) = Cell::make_arrow({Direction::South});
    at(board, 0, 19) = Cell::make_arrow({Direction::South, Direction::West});

    State board_state{board};
    print_board(board_state);
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
