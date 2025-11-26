We need to cut either 2 or 4 squares
- For 2, the corners must be adjacent
- For 4, must include all corner cells

There is a "base" shape with dimensions bounded by the inner corners
of the cut shapes (4) or the inner corners + edge cells in the same 
dimension and width as corners

Valid boxes must have:
- A valid base
- Pair of rectangles adjacent to the base that span to opposite edges must either 
1. Have a combined to-edge width at least the width of base + 2 
-> In this case, the combined width can only be greater in increments of
two cells
2. Have the same to-edge width equal to (combined width of the other pair of rectangles - base width in that dimension) / 2

^ Ok intuition, but we don't need to cut away only rectangles: it's any groupof orthogonal cells with at least one edge cell

Six faces:
Base - 1 
North - 2 
East - 3 
South - 4 
West - 5 
Top - 6 

Can think of this as a graph: if we assign a cell to a face, different movements will either be part of the same face, or move to a different face:
- Assume the base rectangle is static, but we can remove cells from it so long as we add them back from a cell that can reach it
1: N -> 2, E -> 3, S -> 4, W -> 5
2. N -> 2 | 6, E -> 2 (doesn't exceed base cell span E side) | 3 (does), S -> 1 | 2, W -> 2 (doesn't exceed base W span) | 5 (does)
3. N -> 2 | 3,
4. N -> 1 | 4,
5. N -> 2 (exceeds base N span while connected to base) | 

Holdon: think about this in terms of face pairs: think about an arrow on a face pointing away from an adjacent face, four potential per face, i.e. orientation
- More intuitive: can either walk forward or turn, for simplicity assume we turn as if we are inside the box

From the base, we can determine the bounds such that a walk forward in 4 directions will move to another face

From a face adjacent to the base, walks forward may stay in the same face or move to another. The max amount of walks forward in a single face must be equal on opposite faces when extending from the base

If the grid edge is reached and we still want to walk forward on the same face, look for alternative routes

The numbers give us starting face cells: within one king-move distance (including the number cell)
- We can logic our way through some: the 9 is a good starting point: if a number is maxed out, invalidate the remaining adjacent cells
The arrows give non-viable cells
- If two arrows point towards each other and leave one gap, the gap is a face cell
- If a cell has a multi-arrow, the closest cells must have equal distance
Grey squares (7) may be part of the same face group (and thus walkable to each other)
- Squares are walkable to each other if a walk to the cell occupies the adjacent position: think about a cell "wrapping around" to reach the face vs. moving from one face to another
Grey circles may be the same position on opposite face pairs
-> Best to mark these cells first, then fill in the gaps
-> We know we have three rectangle pairs with the same length and width, and number of cells
-> Three dimensional one-hot vector walks 
