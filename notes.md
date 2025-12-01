Grokking

Applying number constraints:
We iterate over r and c -- 0 row is top, 0 col is left: need to adjust my board setup accordingly

Each Cell struct holds a CellContent (an alias for a std::variant object owned by the struct) initialized with an EmptyCell. We have bool helper functions to return whether the variant object holds each of the specific alternatives
-> std::variant is just a class template that represents a type-safe union
-> std::get_if<...>(variant ptr) returns nullptr if the pointer passed as an argument does not point to a type in the variant specified by the template argument

If the board notes have been updated, we return true, else false.
s is a parameter holding the State reference, and s.board gives the std::array holding cell objects
-> Must continue if not a number cell, because these are the only cells we are considering for this function.

When checking connectivity and holes, we use a visited set and a queue. All we have to do is choose a starting cell, add it to visited and queue, then while the queue is not empty, get the queue front Pos, pop from q, inc a visited counter, then for all orthogonal cells, continue if OOB, already in visited, or if not in box, then set the vis for position as true and push the Pos. If after doing this, the visited count is not the same as the entire yes count, we know the box we have isn't fully connected

For holes, we use border cells of the bounding rectangle as potential start points, and explore similarly. If a cell is empty but isn't in the visited empty, we found a hole

So if either check constraints functions return true, we keep checking constraints. If both return false, check connectivity and holes must return true
