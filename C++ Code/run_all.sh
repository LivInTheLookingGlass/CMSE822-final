cd "2D Viscous Burgers Equation Parallel";
python runner.py -O 2 -N_min 5 -N_max 12 -n_min 10000 -n_max 10000 -t_max 8;
python runner.py -O 2 -N_min 5 -N_max 12 -n_min 10000 -n_max 10000 -t_max 8 -u;

cd "../2D Advection Parallel/";
python "../2D Viscous Burger Equation Parallel/runner.py" -O 2 -N_min 5 -N_max 12 -n_min 10000 -n_max 10000 -t_max 8 -f "2D_Advection_Arr_Parallel.cpp";
python "../2D Viscous Burger Equation Parallel/runner.py" -O 2 -N_min 5 -N_max 12 -n_min 10000 -n_max 10000 -t_max 8 -u -f "2D_Advection_Arr_Parallel.cpp";

cd "../2D Burgers Equation Parallel/";
python "../2D Viscous Burgers Equation Parallel/runner.py" -O 2 -N_min 5 -N_max 12 -n_min 10000 -n_max 10000 -t_max 8 -f "2D_Burgers_Arr_Parallel.cpp";
python "../2D Viscous Burgers Equation Parallel/runner.py" -O 2 -N_min 5 -N_max 12 -n_min 10000 -n_max 10000 -t_max 8 -u -f "2D_Burgers_Arr_Parallel.cpp";
