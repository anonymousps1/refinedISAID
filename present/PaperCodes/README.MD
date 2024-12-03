The version of Gurobi we utilized is "Version 10.0.1".
The fold "modelp" is  the MILP model  for PRESENT.
Given a r-round PRESENT, Our modeling sequence is as follows: S--P--K--S--P--k ······ S--P--K--S 


The code to take a backtracking strategy and a restarting strategy is in the fold"./searchcode_backtrack_restart".

Type "make" in the console to generate executable program "genm.run"

If all middle-round key monomial have been obtained , run "python genallkey_12.py" to generate 4096 key monomials for all rounds.

run "./genm.run -r 12" to generate 12-round integral-resistance matrix. 


run "sage determine_rank.sage 12" to solve the rank of 12-round integral-resistance matrix.


