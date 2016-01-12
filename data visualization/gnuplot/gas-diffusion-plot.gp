set term png truecolor enhanced size 640,480
set encoding utf8
set autoscale
set output 'no_no_shift.png'
set title 'NO profile, no shift'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_no_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_no_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_no_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_no_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'o2_no_shift.png'
set title 'O_2 profile, no shift'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_no_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_no_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_no_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_no_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'no_right.png'
set title 'NO profile, right shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_right_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'o2_right.png'
set title 'O_2 profile, right shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_right_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'no_right_down.png'
set title 'NO profile, right shift 20%, down shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_right_down_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_down_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_down_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_down_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'o2_right_down.png'
set title 'O_2 profile, right shift 20%, down shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_right_down_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_right_down_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_right_down_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_right_down_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'no_deform_right.png'
set title 'NO profile, deform 10%, right shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_deformed_right_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'o2_deform_right.png'
set title 'O_2 profile, deform 10%, right shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_deformed_right_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'no_deform_right_down.png'
set title 'NO profile, deform 10%, right shift 20%, down shift 20%'
set ylabel 'NO concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_deformed_right_shift.dat' using 1:2 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:3 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:4 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:5 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'

set autoscale
set output 'o2_deform_right_down.png'
set title 'O_2 profile, deform 10%, right shift 20%, down shift 20%'
set ylabel 'O_2 concentration'
set xlabel 'Radial distance away from vessel center'
plot 'data_small_deformed_right_shift.dat' using 1:6 title '{/Symbol q} = 0' with l lw 2 lt rgb 'red',\
     'data_small_deformed_right_shift.dat' using 1:7 title '{/Symbol q} = 0.5{/Symbol p}' with l lw 2 lt rgb 'green',\
     'data_small_deformed_right_shift.dat' using 1:8 title '{/Symbol q} = {/Symbol p}' with l lw 2 lt rgb 'blue',\
     'data_small_deformed_right_shift.dat' using 1:9 title '{/Symbol q} = 1.5{/Symbol p}' with l lw 2 lt rgb 'cyan'
