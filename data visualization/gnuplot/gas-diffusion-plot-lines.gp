reset
set term png truecolor enhanced size 1200,480
set encoding utf8
set autoscale
set output 'lines_no_no_shift.png'
set title 'NO profile, no shift'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
set arrow 1 from 25,graph(0,0) to 25,graph(1,1) nohead lt rgb 'black'
set arrow 2 from 20.7,graph(0,0) to 20.7,graph(1,1) nohead lt rgb 'red'
set arrow 3 from 20.7,graph(0,0) to 20.7,graph(1,1) nohead lt rgb 'green'
set arrow 4 from 20.7,graph(0,0) to 20.7,graph(1,1) nohead lt rgb 'blue'
set arrow 5 from 20.7,graph(0,0) to 20.7,graph(1,1) nohead lt rgb 'cyan'
plot 'data_small_no_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_no_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_no_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_no_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_o2_no_shift.png'
set title 'O_2 profile, no shift'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_no_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_no_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_no_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_no_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_no_right.png'
set title 'NO profile, right shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
set arrow 2 from 21.5,graph(0,0) to 21.5,graph(1,1) nohead lt rgb 'red'
set arrow 3 from 20.6,graph(0,0) to 20.6,graph(1,1) nohead lt rgb 'green'
set arrow 4 from 19.8,graph(0,0) to 19.8,graph(1,1) nohead lt rgb 'blue'
set arrow 5 from 20.6,graph(0,0) to 20.6,graph(1,1) nohead lt rgb 'cyan'
plot 'data_small_right_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_o2_right.png'
set title 'O_2 profile, right shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_right_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_no_right_down.png'
set title 'NO profile, right shift 20%, down shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
set arrow 2 from 21.5,graph(0,0) to 21.5,graph(1,1) nohead lt rgb 'red'
set arrow 3 from 19.8,graph(0,0) to 19.8,graph(1,1) nohead lt rgb 'green'
set arrow 4 from 19.8,graph(0,0) to 19.8,graph(1,1) nohead lt rgb 'blue'
set arrow 5 from 21.5,graph(0,0) to 21.5,graph(1,1) nohead lt rgb 'cyan'
plot 'data_small_right_down_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_down_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_down_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_down_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_o2_right_down.png'
set title 'O_2 profile, right shift 20%, down shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_right_down_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_down_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_down_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_down_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_no_deform_right.png'
set title 'NO profile, deform 10%, right shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
set arrow 2 from 21.5,graph(0,0) to 21.5,graph(1,1) nohead lt rgb 'red'
set arrow 3 from 18.6,graph(0,0) to 18.6,graph(1,1) nohead lt rgb 'green'
set arrow 4 from 19.8,graph(0,0) to 19.8,graph(1,1) nohead lt rgb 'blue'
set arrow 5 from 22.7,graph(0,0) to 22.7,graph(1,1) nohead lt rgb 'cyan'
plot 'data_small_deformed_right_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_o2_deform_right.png'
set title 'O_2 profile, deform 10%, right shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_deformed_right_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_no_deform_right_down.png'
set title 'NO profile, deform 10%, right shift 20%, down shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
set arrow 2 from 21.5,graph(0,0) to 21.5,graph(1,1) nohead lt rgb 'red'
set arrow 3 from 17.8,graph(0,0) to 17.8,graph(1,1) nohead lt rgb 'green'
set arrow 4 from 19.8,graph(0,0) to 19.8,graph(1,1) nohead lt rgb 'blue'
set arrow 5 from 23.7,graph(0,0) to 23.7,graph(1,1) nohead lt rgb 'cyan'
plot 'data_small_deformed_right_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'lines_o2_deform_right_down.png'
set title 'O_2 profile, deform 10%, right shift 20%, down shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_deformed_right_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'
