set terminal  postscript eps size 3.5,2.5 enhanced color  linewidth 2  font 'Helvetica,12'

#comentarios "dentro" do gnuplot iniciam com #

 set output 'sigtotal.eps' 
#Tamanho do gráfico (cortamos o excesso com o comando crop contigo em run0.sh)   
  set size   0.95,0.9
#Divisao de n partes nos eixos x e y    
  set mxtics 4
  set mytics 5

#Intervalo de x e y , xtics=1 -> 1,2,3 , xtics=0.5 ->1.0,1.5,2.0	 
   #set xtics  1
   #set ytics 0.1

# Ativando a rede
#  set grid
#Intervalo do gráfico    
 # set yrange [0.01:1.0]
 #  set xrange [-0.8:0.8]
# Escala logarítmica em x e y.  
#   set logscale x
 #  set logscale y

#  set format y '10^{%T}'
#  set format x '10^{%T}'
 #set label 1  'Data from Tasso Experiment' at 0.4,0.002  center
 set xlabel " {/Symbol \326}s  [GeV] " font "Times-Italic,10" 
 set ylabel "{/Symbol s} [nb]" font "Times-Italic,10"
  set title ' Bhabha elastic scattering ' font "Times-Bold,14" 
  set key top right
 # set fit quiet
 # f(x)=a*x**(-b)
 # fit f(x) 'WFIT.dat' via a,b
 # set print 'fitp.dat'
 # print b 
 # set print
  plot 'BHABHA_TOT.dat' using 1:2 with lines lc rgb 'black' linetype 5 linewidth 1 title 'QED Theory',\
       'sigtotalexp.dat'	using 1:2:3 with yerrorbars  lc rgb 'red' linetype 1 title'Tasso Experiment'
	
	
