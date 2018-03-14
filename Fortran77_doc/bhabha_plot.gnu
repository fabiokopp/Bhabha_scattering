set terminal  postscript eps size 3.5,2.5 enhanced color  linewidth 2  font 'Helvetica,12'

#comentarios "dentro" do gnuplot iniciam com #

 set output 'bhabha.eps' 
#Tamanho do gráfico (cortamos o excesso com o comando crop contigo em run0.sh)   
  set size   0.95,0.9
#Divisao de n partes nos eixos x e y    
  #set mxtics 2
  set mytics 5

#Intervalo de x e y , xtics=1 -> 1,2,3 , xtics=0.5 ->1.0,1.5,2.0	 
   set xtics  0.1
   #set ytics 0.1

# Ativando a rede
#  set grid
#Intervalo do gráfico    
 # set yrange [0.01:1.0]
   set xrange [-0.8:0.8]
# Escala logarítmica em x e y.  
#   set logscale x
   set logscale y

#  set format y '10^{%T}'
#  set format x '10^{%T}'
 set label 1  'Data from Tasso Experiment' at 0.4,0.002  center
 set xlabel "cos {/Symbol q} " font "Times-Italic,10" 
 set ylabel "d{/Symbol s}/d{/Symbol W} [nb/sterad]" font "Times-Italic,10"
  set title ' Bhabha elastic scattering ' font "Times-Bold,14" 
  set key top left
 # set fit quiet
 # f(x)=a*x**(-b)
 # fit f(x) 'WFIT.dat' via a,b
 # set print 'fitp.dat'
 # print b 
 # set print
 p1=14.0**2.0
 p2=22.0**2.0
 p3=34.8**2.0
 p4=38.3**2.0
 p5=43.6**2.0
  plot 'BHABHA.dat' using 1:2 with lines lc rgb 'black' linetype 5 linewidth 1 title 'QED Theory',\
	'experimental-data.txt'	using 1:($4/p1):($5/p1) with yerrorbars  lc rgb 'red' linetype 1 title'{/Symbol \326}s=14.0 GeV',\
	'experimental-data.txt'	using 1:($7/p2):($8/p2) with yerrorbars  lc rgb 'blue' linetype 8 title'{/Symbol \326}s=22.0 GeV',\
       	'experimental-data.txt'	using 1:($10/p3):($11/p3) with yerrorbars  lc rgb 'cyan' linetype 7 title'{/Symbol \326}s=34.8 GeV',\
	'experimental-data.txt'	using 1:($13/p4):($14/p4) with yerrorbars  lc rgb 'brown' linetype 4 title'{/Symbol \326}s=38.3 GeV',\
	'experimental-data.txt'	using 1:($16/p5):($17/p5) with yerrorbars  lc rgb 'yellow' linetype 5 title'{/Symbol \326}s=43.6 GeV'
   
 
	
