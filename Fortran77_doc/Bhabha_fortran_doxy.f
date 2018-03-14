!>       Programa do processo Bhabha. Calcula a seção de choque
!!       diferencial para uma determinada energia de centro de massa, bem 
!!       como a seção de choque total em termos da energia de centro de massa. \n
!!	A compilação deve ser feita segundo o comando: gfortran -o Bhabha.exe Bhabha_fortran_doxy.f subroutines.f \n
!!	Ou a versão sem documentação incluída e com as rotinas no mesmo arquivo: gfortran -o Bhabha.exe Bhabha_fortran_sem_doxy.f \n 
!!	A execução ocorre via: ./Bhabha.exe \n
!!	Os gráficos dependem do programa gnuplot (cujos arquivos possuem extensão .gnu) para serem gerados. Em distribuições Debian, a instalação ocorre via:\n
!!	sudo apt-get install gnuplot \n
!!	Autor: \f[ \bf{Fabio \  Kopp  } , \ Instituto \ de  \ Física \ (UFRGS), \ RS , \ Brasil. \f] \n
!!	Email: \f[ \bf{  \ fabio.kopp@ufrgs.br} \f]	\n
!!	Versões: 0.1a  - (13/03/2018).  \n
	PROGRAM BHABHA_SCATTERING
	IMPLICIT NONE
	DOUBLE PRECISION DSIGMA,Ecm,X,Ecm0,SIGTOTAL,S2
	DOUBLE PRECISION chi, chired,S20
	INTEGER I

c      Cria dois arquivos .dat utilizados para salvar os resultados de dsigma/domega e sigma_total.
	OPEN(UNIT=10, FILE="BHABHA.dat", STATUS="UNKNOWN")
	OPEN(UNIT=20, FILE="BHABHA_TOT.dat", STATUS="UNKNOWN")

c      Energia de centro de massa em GeV
	WRITE(*,*) "Entre com Sqrt(S) [GeV]"
	READ(*,*) S2

c      x = cos(theta)	 
	DO I=-80,80,1
	X=I*0.01D0
c      Imprimi na tela o cos(theta) e dsigma/domega	
	PRINT*, "cos(theta) ,  dsigma/domega (nb/sterad)"
	PRINT*, X, DSIGMA(S2,X)

c      Salva o cos(theta) e dsigma/domega no arquivo BHABHA.dat	
	WRITE(10,*) X,DSIGMA(S2,X) 	
	END DO

c      Fecha o arquivo Bhabha.dat file
	CLOSE(10)

	DO I=10,50,1
	S20=I*1.0D0
	PRINT*, "S^(1/2)(GeV)  , SIGMATOTAL (nb)"
	PRINT*, S2, SIGTOTAL(S2)
c      Salva a S2 e SIGTOTAL no arquivo BHABHA-TOT.dat	
	WRITE(20,*)  S20, SIGTOTAL(S20)
	END DO
	CLOSE(20)

c	Chama a subrotina a qual faz a analise estatistica
	call chisquared(19,0,chi,chired,34.8d0)
	
c      Executa no terminal o gnuplot e abre os graficos com o programa evince 
	CALL SYSTEM("gnuplot bhabha_plot.gnu& evince bhabha.eps&")
	CALL SYSTEM("gnuplot sigtotal.gnu& evince sigtotal.eps&")

	
	END PROGRAM BHABHA_SCATTERING	
	
!>      Sessão de choque diferencial  \f[ \frac{d\sigma}{d\Omega}(S2= \sqrt(s),x = cos(\theta))[\frac{nb}{sterad}] \f]
	DOUBLE PRECISION FUNCTION DSIGMA(S2,X)
	DOUBLE PRECISION ALPHAQED2,X,DS,FACTOR1
	DOUBLE PRECISION CF,S,S2
	DOUBLE PRECISION FACTOR2, FACTOR3,Ecm
	INTEGER I
	
c      Constante de estrutura fina ao quadrado (PDG(2016))	
	ALPHAQED2=(1.0D0/137.035999074d0)**2.0

C 	DEFINITION OF CHANNEL S
C	S=4.0*(Ecm**2)  -> Ecm**2 = S/4


C 	CONVERSION FACTOR FROM GeV^-2 TO NB 
C       CF=0.389*d6 
	S=(S2)**2.0
	CF=0.389d6
	FACTOR1=((1.0D0 + X)**2.0 + 4.0D0)/
     &((1.0D0-X)**2.0)
	FACTOR2=(1.0D0+X**2.0)/2.0D0
	FACTOR3=-((1.0D0+X)**2.0)/(1.0D0-X)	
	DSIGMA=(ALPHAQED2/(2.0*(S)))*(FACTOR1+FACTOR2+
     & FACTOR3)*CF
	RETURN 
	END

!>      Função auxiliar para a integração da seção de choque diferencial no ângulo sólido
	DOUBLE PRECISION FUNCTION  FAUX(X)
	IMPLICIT NONE
	DOUBLE PRECISION X,S2,DSIGMA,PI,T1
	COMMON/EN/S2
	PI=3.141592653d0
!>      \f[ \frac{d\sigma}{d\Omega} \cdot 2.0 \ \pi \ \cdot sin(\theta) \f]
c       Integration over the solid angle, e.g 2*PI*Sin(theta)*(dsigma/domega) d(theta)
	T1=(2.0*PI*(sqrt(1.0-X**2)))
	FAUX=DSIGMA(S2,X)*T1
	RETURN
	END
	
!>      Seção de choque total \f[ \sigma( \ S20 = \sqrt(s) \ ) [nb] \f], resultado da integração da função FAUX(X) usando
!! a rotina de integração Simpson(função a ser integrada, limite inferior, limite superior, número de intervaloes)	
	DOUBLE PRECISION FUNCTION SIGTOTAL(S20)
	IMPLICIT NONE
	EXTERNAL f, DSIGMA,FAUX
	DOUBLE PRECISION S2,simpson,FAUX,S20
	COMMON/EN/S2
	S2=S20
	
!>	\f[ \sigma[nb] = \int_{-0.84}^{0.84}\! \frac{d\sigma}{d\Omega} \cdot 2.0 \ \pi \ \cdot sin(\theta) \, \mathrm{d}\theta \f]
 
c       If we consider a deviation in the theta about 5%, then 0.8->0.84 and the plot matches with the experimental plot 
c presented in the paper.
c       However, is we take the exact interval integration, the curve moves to below of experimental points. 
c                simpson(function to be integrated,inferior limit, superior limit, number of intervals)
	SIGTOTAL=simpson(FAUX,-0.84d0,0.84d0,100)
	RETURN
	END

c       N is the number of points (data file),v is the number of variables in the model, chi is the chi-squared value , chired is the chi-squared reduced and si2 is energy of center of mass. 
	subroutine chisquared(N,v,chi,chired,si2)
	implicit none
	integer N,I,v
	double precision chi, chired,x(N),y(N),d2y(N)
	double precision sol,si2
	double precision yth(N),s,sf(N)
	double precision C(18,N),DSIGMA
!>      Entrada: \n
!!	\n
!!      N é o número de pontos no arquivo de dados. \n
!! 	v é o número de variáveis do modelo. \n
!!	si2 é a energia de centro de massa em GeV.	\n
!!      \n
!!	Saída: \n 	
!!	chi é o valor do Qui-quadrado em uma dada energia de centro de massa. \n
!!      chired é o Qui-quadrado reduzidos. \n
!!	\n
!!      Para calcular o Qui-quadrado de uma energia de centro de massa diferente de 34.8 GeV,
!!	altere as colunas de y(N)=C(?,N) e d2y(N)=sqrt((C(?,N)/s)**2) para os valores experimentais 
!!      relacioados a outra energia de centro de massa. Lembrando que y(N) é o valor da seção de choque diferencial e d2y o erro 
!! 	desta somado em quadratura.

	open(UNIT=30, FILE="experimental_data.dat", STATUS="UNKNOWN")

	s=si2**2.0
	sol=0.0d0

	do I=1,N,1
	read(30,*) C(1,N),C(2,N),C(3,N),C(4,N),C(5,N),C(6,N), 
     &C(7,N),C(8,N),C(9,N),C(10,N),C(11,N),C(12,N),C(13,N),
     &C(14,N),C(15,N),C(16,N),C(17,N) , C(18,N)

c       y and d2y is referent to a s^1/2=34.8 GeV	
c 	d2y is calculated in quadrature ( we  add the two erros in quadrature)
	x(N)=C(1,N)  
        y(N)=C(10,N)/s
	d2y(N)=sqrt((C(11,N)/s)**2)
c      yth is the dsigma/domega in which we evaluated at experimental values of cos(theta)
	yth(N)=DSIGMA(si2,C(1,N))
	sol=sol+ ((y(N)-yth(N))/d2y(N))**2
c	write(*,*), x(N) , y(N) , d2y(N), yth(N) 
	end do
	close(30)

	chi=sol 
	chired=sol/(N-v)	
	

	write(*,*)"--------------statistical analysis--------------"
	write(*,'(1x,a,4x,f5.2,2x,a)'), "s^1/2 = " , si2 , "GeV"
	write(*,'(1x,a,1x,f15.7)') "Chi-Square =" , chi 
	write(*,'(1x,a,1x,i3)') "Ndf = " , N	
	write(*,'(1x,a,1x,i3)') "Number of variables to be adjusted = " ,
     & v	
	write(*,'(1x,a,1x,f15.7)') "Chi-square/ndf [reduced] =" ,
     & chired 
	call probchi2(chi,(N-v)*1.0d0)
	write(*,*)"-----------------------------------------------"
	return 
	end
