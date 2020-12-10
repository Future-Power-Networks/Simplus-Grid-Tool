a = sym('a');
b = sym('b');
c = sym('c');

u = [a; b; c];
% u = [1;2;3];
t1 = SimplexPS.abc2alphabeta(u,'ZeroSequence',1)
t2 = SimplexPS.alphabeta2abc(t1)
t3 = simplify(t2)