
############################
#
#
#   Program : simulations associated to the paper arXiv:2203.16160v1 
#  
#     On the long time behaviour of single stochastic Hodgkin-Huxley neurons 
#     with constant signal, and a construction 
#     of circuits of interacting neurons 
#     showing self-organized rhythmic oscillations 
#
#   submitted to Mathematical Neurosciences and Applications (MNA) 
#
#
#   Reinhard Hoepfner, Institute of Mathematics, University of Mainz 
#
#
#   Update 26.08.22  
#
#
#
#   0) : Preliminaries 
#
#   I) : Single stochastic Hodgkin Huxley neuron with constant signal 
# 
#   II) : Transfer functions   
#   
#   III) : Circuits  
#
#
############################






#############################
#############################
#
#   0) : preliminaries 
#
############################
############################



############################
#
#   0 a) step size and time horizon for euler schemes 
#
###########################
 
 
delta <- 0.001 ; 
# delta <- 0.004 ; 
  # step size for euler schemata 
  # for simulation of stochastic processes  

T_lim <- 600 ; 
# T_lim <- 1500 ; 
  # time horizon for euler schemata 
  # for simulation of stochastic processes  

tt <- seq( 0, T_lim, delta ) ; 
length(tt) ; 
  # the time grid 



#############################
#
#   0 b) the functions govering the internal variables 
#        in the classical deterministic Hodgkin-Huxley model 
#        (constants as in the book Izhikevich 2007) 
#
###############################


alpha_n <- function(v){ 
   w <- ifelse( v==10, 10.0001, v ) ; # avoid division by 0 
   0.01*(10-w) / (exp(0.1*(10-w))-1) 
   } ; 

beta_n <- function(v){ 0.125*exp(-v/80) } ; 

alpha_m <- function(v){    
   w <- ifelse( v==25, 25.0001, v ) ; # avoid division by 0 
   0.1*(25-w) / (exp(0.1*(25-w))-1) 
   } ; 

beta_m <- function(v){ 4*exp( -v/18 ) } ; 

alpha_h <- function(v){ 0.07*exp( -v/20 ) } ;        

beta_h <- function(v){ 1 / ( exp(0.1*(30-v)) + 1 ) } ;  




###############################
#
#   0 c)  increments in the deterministic Hodgkin Huxley model when step size is  delta  
#
###############################


   #  membrane potential increment 
   #  due to the function  F  in the HH model (constants as in the book Izhikevich 2007) 
   #  when step size is  delta  

d_IF <- function( v, n, m, h, delta ) { 
  ( 36 * n^4 * (v+12) + 120 * m^3*h * (v-120) + 0.3 * (v-10.6) ) * delta 
  } ; # function calculating increments  F(V_t,n_t,m_t,h_t) dt 


   # increments for the internal variables  m , n , h 
   #  when step size is  delta 

d_n  <- function( v, n, delta ){ ( alpha_n(v) * (1-n) - beta_n(v) * n ) * delta } ; 

d_m  <- function( v, m, delta ){ ( alpha_m(v) * (1-m) - beta_m(v) * m ) * delta } ;

d_h  <- function( v, h, delta ){ (alpha_h(v) * (1-h) - beta_h(v) * h ) * delta } ;


   
   
##############################
#
#   0 d) increments of the Ornstein Uhlenbeck process 
#        with backdriving force  tau > 0  and volatility  sigma > 0  
#        (in absence of any signal, i.e. mean-reverting to  0 ) 
#
#        OU invariant law is centred normal with variance  gamma^2/2  ,  
#        gamma := sigma / sqrt(tau)  
#
##############################


   # chose backdriving force  tau  and volatility  sigma  
sigma <- 1.5 ; tau <- 1.4 ; 
# sigma <- 1 ; tau <- 0.7 ; 

   # then define  gamma   by  
gamma <- sigma / sqrt(tau) ; 
  
  
   # increments in euler schemes for the Ornstein Uhlenbeck process   phi  
   # (without signal) when step size is  delta  
   
d_phi <- function( phi_alt, delta, gamma, tau ){  
   -tau * phi_alt * delta + gamma * sqrt(tau*delta) * rnorm(1,0,1)   
   } ; # function calculating increments in euler schemes 
   #  when step size is  delta  
 
 
 
 
 
 
#####################################
#####################################
#
#   I) : Single stochastic Hodgkin Huxley neuron with constant signal 
#
#   (this corresponds to section 2 
#   in the paper arXiv:2203.16160v1 submitted to MNA)
#
#####################################
#####################################


   # in view of later use with time-varying input, 
   # we consider a single stochastic Hodgkin Huxley neuron 
   # where the first line of the Hodgkin Huxley equations is written as  
   #
   # $ dV_t = N0_input_faktor dt + dphi_t - F(V_t, n_t, n_t, m_t, h_t ) dt $    
   # 
   # and simulate the following variables of the stochastic Hodgkin Huxley neuron 
   # as considered in section 2 of the paper arXiv:2203.16160v1 submitted to MNA :    
   #
   #  a) trajectories of the biological variables : 
   #     membrane potential t  -->  V_t 
   #     and gating variables   t  -->  n_t,  m_t,  h_t    
   #
   #  b) Ornstein Uhlenbeck trajectory   t  -->  phi_t   
   #
   #  c) the spiking pattern  spike_id  , 
   #     i.e. the point process of spike times 
   #     where we define the beginning of a spike 
   #     by upcrossing of the  m - variable over the  h -variable 
   #
   #  d) trajectory  t  -->  U_t   of the output of this neuron, 
   #     defined as in section 2.3 of the paper arXiv:2203.16160v1   
   #     as solution to 
   #
   #     $ dU_t = -\lambda U_t + c dN_t $ 
   #
   #     where  N  is the counting process associated to the point process of spike times 
   #     (in section 2.3 of the paper arXiv:2203.16160v1 submitted to MNA, we put  c := 1 , 
   #     and write  c_1  instead of  lambda )  
   #
   #  the eight variables specified in a) -- d)  will be calculated below 
   #  as a function of time, using an euler scheme with step size  delta    
    


#################################
#
#   I i ) : parameters for the output process 
#
#################################


lambda <- 0.02 ; 
   # decay parameter for the output process 
   # in section 2.3 of the paper arXiv:2203.16160v1 submitted to MNA, 
   # we write  c_1  instead of  lambda    
   
   # caveat : recalibration of the decay parameter 
   # according to 5.3 in the paper arXiv:2203.16160v1 
   # might be necessary once we know 
   # a median of interspike times under regular spiking  
   # as indicated in  I vi)  below 
   
   
c <- 1 ; 
   # in section 2.3 of the paper arXiv:2203.16160v1  
   # we always consider  c := 1   
   
   
 
#################################
#
#   I ii ) : choice of values for the signal 
#
#################################


  # in equations (8)-(11) of the paper arXiv:2203.16160v1 submitted to MNA, 
  # values for a time-constant signal in a stochastic Hodgkin-Huxley system 
  # are denoted by  theta ; in view of later use with time-varying input,  
  # we write   N0_input_faktor  instead of   theta  
  
   
N0_input_faktor <- 10 ; 
  # a deterministic neuron with signal  theta = 10  
  # would be attracted to a stable orbit and thus  
  # regularly spiking 

# N0_input_faktor <- 4 ; 
  # a deterministic neuron with signal  theta = 4  
  # would be attracted to the stable equilibrium , 
  # i.e. would eventually be quiet 

# N0_input_faktor <- 6 ; 
  # for a deterministic neuron, signal  theta = 6   
  # would belong to the bistability interval 
  # described in section 1 of the paper arXiv:2203.16160v1 : 
  # depending on the initial conditions, 
  # trajectories would be attracted 
  # either by the stable orbit or by the stable equilibrium  

  # in stochastic Hodgkin Huxley systems with signal  theta = 6 , 
  # in the long run, 
  # long periods of seemingly regular spiking 
  # should alternate with long quiet periods 



#############################################
#
#  I iii)  preparing vectors in order to store 
#          the eight variables specified in a) -- d) as a function of time,  
#          random choice of starting values  
#
#############################################


N0_input <- N0_input_faktor * rep( 1, length(tt) ) ; 
   # signal constant in time 

N0_phi <- rep(0,length(tt)) ;
N0_phi[1] <- rnorm( 1, 0, (gamma/sqrt(2)) )  ; 
   # for Ornstein Uhlenbeck process  phi  running stationary  

N0_v <- rep(0,length(tt)) ; 
N0_v[1] <- runif( 1, -12, 120 )  ; 
   # for membrane potential  V 

N0_n <- rep(0,length(tt)) ; 
N0_n[1] <- runif( 1, 0, 1 ) ; 
   # for gating variable  n  

N0_m <- rep(0,length(tt)) ; 
N0_m[1] <- runif( 1, 0, 1 ) ; 
   # for gating variable  m  

N0_h <- rep(0,length(tt)) ; 
N0_h[1] <- runif( 1, 0, 1 ) ;
   # for gating variable  h  

N0_spike <- rep(0,length(tt)) ; 
   # for spiking pattern 

N0_output <- rep( 0, length(tt) ) ; 
  # for output  U  



#############################################
#
#  I iv)  simulation 
#
#         single stochastic Hodgkin Huxley neuron with constant signal   
#         the first line of the Hodgkin Huxley equations is 
#
#         $ dV_t = N0_input_faktor dt + dphi_t - F(V_t, n_t, n_t, m_t, h_t ) dt $  
#
#         via euler scheme of step size  delta  
#
#         some related graphics 
#
#############################################
 
   
j <- 2 ; while( j <= length(tt) ){
   # 
   input_alt <- N0_input[j-1] ; 
   phi_alt <- N0_phi[j-1] ; 
   v_alt <- N0_v[j-1] ; 
   n_alt <- N0_n[j-1] ; 
   m_alt <- N0_m[j-1] ; 
   h_alt <- N0_h[j-1] ; 
   output_alt <- N0_output[j-1] ; 
   # 
   dphi <- d_phi( phi_alt, delta, gamma, tau ) ; 
   dif <- d_IF( v_alt, n_alt, m_alt, h_alt, delta ) ; 
   v_neu <- v_alt + ( input_alt*delta + dphi ) - dif ; 
   n_neu <- n_alt + d_n( v_alt, n_alt, delta ) ; 
   m_neu <- m_alt + d_m( v_alt, m_alt, delta ) ; 
   h_neu <- h_alt + d_h( v_alt, h_alt, delta ) ; 
   dN <- ifelse( m_alt<h_alt && m_neu>h_neu, 1, 0 ) ; 
   output_neu <-  output_alt - lambda*output_alt*delta + c*dN    ; 
   #
   #N0_input[j] <- ..... 
   N0_phi[j] <- phi_alt + dphi ; 
   N0_v[j] <- v_neu ; 
   N0_n[j] <- n_neu ; 
   N0_m[j] <- m_neu ; 
   N0_h[j] <- h_neu ; 
   N0_spike[j] <- dN ; 
   N0_output[j] <-  output_neu ; 
   #cat( j, "\n" ) ; 
   j <- j+1 ; 
   # 
   } ; # ende while 
   # this loop calculates the eight variables in  a) -- d)  as a function of time 




###############################################################
#
#   I v) test graphics showing the eight variables 
#        caclulated in  a) -- d)  above 
#        (single stochastic Hodgkin Huxley neuron with constant signal) 
#
###############################################################
   
   
leg <- "membrane potential  V  when  dY_t = !!!! dt + dX_t " ; 
leg1 <- "OU process  X  ( with  tau = §§§§ , sigma = %%%% , gamma = &&&& )" ; 
leg1 <- gsub( "§§§§", round(tau,2), leg1 ) ; 
leg1 <- gsub( "&&&&", round(gamma,2), leg1 ) ; 
leg1 <- gsub( "%%%%", round(sigma,2), leg1 ) ; 
leg3 <- " output process  U " ; 
leg2 <- "gating variables  n (green) , m (blue), h (magenta)" ; 
leg <- gsub( "!!!!", N0_input[1], leg ) ; leg ; 
q1 <- qnorm( 0.01, 0, gamma/sqrt(2) ) ; # invariant law of OU process, a lower 1% quantile  
q2 <- qnorm( 0.99, 0, gamma/sqrt(2) ) ; # invariant law of OU process, an upper 1% quantile  
#
par(mfrow=c(2,2) ) ; 
#
plot( range(tt), c(-20,120), type="n", main=leg, ylab="", xlab="time" ) ; 
for( j in 1:length(tt) ) { if ( N0_spike[j] == 1 ) abline( v = tt[j], lty=2 ) ; } ; # spiking times 
lines( tt, N0_v, col=2 ) ; # membrane potential  V 
#
plot( range(tt), range(N0_phi), type="n", main=leg1, ylab="", xlab="time" ) ;
lines( tt, N0_phi ) ; 
abline( h=q1, lty=2 ) ; 
abline( h=q2, lty=2 ) ; # OU trajectory with quantiles of invariant law 
#
xleg <- "blue : m , green : n , magenta : h"
plot( range(tt), range(-0.1,1.1), type="n", main=leg2, ylab="", xlab="" ) ; 
for( j in 1:length(tt) ) { if ( N0_spike[j] == 1 ) abline( v = tt[j], lty=2 ) ; } ; # spiking times 
lines( tt, N0_n, col=3 ) ; # gating variable  n  
lines( tt, N0_m, col=4 ) ; # gating variable  m 
lines( tt, N0_h, col=6 ) ; # gating variable  h 
#
plot( range(tt), range(N0_output), type="n", main=leg3, xlab="", ylab="" ) ; 
lines( tt, N0_output, col=4 ) ; # output  U  
# 
par(mfrow=c(1,1) ) ;  


   # suggestion : to see what happens when the signal belongs 
   # to the bistability interval of the deterministic HH model, 
   # run the program from the line ' N0_input_faktor <- 6 ; ' 
   # to the present graphics  

   # in this way we have generated the graphics  figures 1, 2, 6  
   # in the paper arXiv:2203.16160v1 submitted to MNA 
   # for the single stochastic Hodgkin Huxley neuron with constant signal 




########################################
#
#   I vi) : preliminaries on spiking patterns and  median of interspike times 
#           when spiking is regular enough 
#
#########################################



   # here we assume that 
   # a sufficient number of sufficiently regularly spaced spikes 
   # could be observed 
   # in the preceding simulation 
   # of a stochastic Hodgkin Huxley neuron :  
   

   # visual check : is there some regularity in the spiking pattern ? 
spikezeiten <- tt[ N0_spike > rep(0,length(N0_spike)) ] ; 
spikezeiten ;   # identification of spike times 
length(spikezeiten) ; # the total number of spikes  
plot( c(0,T_lim), c(-1,1), type="n", main="the spiking pattern", xlab="", ylab="" ) ; 
for( j in 1:length(spikezeiten) ) abline( v=spikezeiten[j], lty=2 ) ; 
  
   # consider the interspike times and their median : 
   
spikezeitendiff <- spikezeiten[-1] - spikezeiten[-length(spikezeiten)] ; 
spikezeitendiff ; 
length(spikezeitendiff) ; 
sort(spikezeitendiff) ; 
Delta <- median(  spikezeitendiff ) ; 
Delta ; # the median  Delta  of the interspike times will play 
   # --under certain assumptions -- an important role in the sequel 
   # cf. formulae (46) and (47) of the paper arXiv:2203.16160v1 submitted to MNA 


   # the vector  spikezeitendiff  of interspike times 
   # is the basis for statistical analysis of interspike times, 
   # e.g. empirical distribution functions, empirical laplace transforms 
   # as in section 3 of the paper arXiv:2203.16160v1  
   
   
   # now consider the median of the interspike times and ratios 
   # ( upper minus lower alpha quantile ) / median 
   # as in section 4 of the paper arXiv:2203.16160v1 
   # in order to judge wether or not 
   # spiking is regular enough in the sense of definition 4.1 of the paper : 
   
alpha1 <- 0.1 ; 
# alpha1 <- 0.05 ; 
# alpha1 <- 0.25 ; 
ratio_alpha1 <- ( quantile( spikezeitendiff, 1-alpha1 ) - quantile( spikezeitendiff, alpha1 ) ) / Delta ; 
ratio_alpha1 ; 

   

   # in case of sufficiently regular spiking : 
   # plot orbits  V  against  n  (for a regularly spiking neuron) 
   # and visualize the definition of spiking times 
   # by upcrossing of  m  over  h 
   # as in formula (16) of the paper arXiv:2203.16160v1 submitted to MNA    
   
   # we show only the second half of the trajectory  
   
hleg <- "orbits  V  against  n  (where we put a dot every 0.1 sec) and  beginning of the spikes as defined in the paper" ; 
xleg <- "membrane potential  V  (parameters of the stochastic Hodgkin Huxley neuron :  theta = §§§§ , sigma = %%%% , tau = &&&& )" ; 
xleg <- gsub( "§§§§", N0_input_faktor, xleg ) ; 
xleg <- gsub( "%%%%", round(sigma,2), xleg ) ; 
xleg <- gsub( "&&&&", round(tau,2), xleg ) ; 
yleg <- "gating variable  n " ; 
plot( c(0.3,0.85), c(-30,120), type="n", main=hleg, xlab=xleg, ylab=yleg ) ; 
i <- round( 0.5*length(tt) ) ; while( i < length(tt) ){ 
   points( N0_n[i], N0_v[i], cex=0.1 ) ;
   i <- i+100 ; 
   } ; # for delta = 0.001 , we put a dot every 0.1 sec 
i <- round( 0.5*length(tt) ) ; while( i < length(tt) ){ 
   if( N0_spike[i] == 1){ points( N0_n[i], N0_v[i], col=2 ) ; } ; 
   i <- i+1 ; 
   } ; # this situates the beginning of the spikes 
   # on the orbit V  against  n 




########################################
#
#   I vii) : visual check : properties of the output process 
#            in the long run 
#            when spikes occur regularly enough 
#
#########################################

   # under the following assumption : 
   
   # inspection of the spiking pattern in  I vi) 
   # did show enough regularity 


   # we visualize the benchmarks proposed 
   # for the range of the output process in the long run 
   # using the median  Delta  of the interspike times 
   # as in formulae (37) and (43) 
   # or in formulae (46) and (47) 
   # of the paper arXiv:2203.16160v1 submitted to MNA 
   
   
u_inf <- exp( - lambda *Delta ) / ( 1 - exp( - lambda *Delta ) ) ; 
u_sup <- u_inf + 1 ; 
u_inf ;  # in the long run : benchmark 
   # for the state of the output process immediately before a spike 
u_sup ; # in the long run : benchmark 
   # for the state of the output process immediately after a spike 


   # compare the trajectory of the output process 
   # to the benchmarks  u_sup and  u_inf  :  
   
output_after_spike <- N0_output[ N0_spike > rep(0,length(N0_spike)) ] ; 
output_after_spike ; 
sort( output_after_spike ) ;  # caveat : there is a starting period 
median( output_after_spike ) ; # compare to  u_sup 
   #
output_before_spike <- output_after_spike[-1]- rep( 1, length(output_after_spike[-1]) ) ; 
output_before_spike ; 
sort( output_before_spike ) ;  # caveat : there is a starting period 
median(  output_before_spike ) ; # compare to   u_inf 
   # 
leg4 <- "output process  U  ( driving OU with parameters tau = §§§§ , sigma = %%%% )  and the benchmarks  u_inf , u_sup " ; 
leg4 <- gsub( "§§§§", round(tau,2), leg4 ) ; 
leg4 <- gsub( "%%%%", round(sigma,2), leg4 ) ; 
plot( range(tt), range(N0_output), type="n", main=leg4, xlab="", ylab="" ) ; 
lines( tt, N0_output, col=4 ) ; 
abline( h = u_inf , col=2, lty=2 ) ; 
abline( h = u_sup , col=2, lty=2 ) ; 



   # CAVEAT : one may have to go back to  I i) 
   # and recalibrate the decay parameter for the output process 
   # of the regularly spiking neuron 
   # according to 5.3 of the paper arXiv:2203.16160v1 
   
   
  


###########################################
###########################################
#
#   II)  Simulation of circuits of interacting Hodgkin Huxley neurons 
#
#        as constructed in section 5 
#        of the paper arXiv:2203.16160v1 submitted to MNA 
#
#        (three blocs of equal length, 
#        interaction clockwise along the circuit, 
#        excitation within blocs, 
#        inhibition from every bloc to its successor) 
#
#        the simulation will be done in II v)  
#
#        preparations in  II i) -- iv) 
#
############################################
############################################


   #   i)    definition of transfer functions 
   #
   #   ii)   chosing the size of the circuit, preparations 
   #         (in the framework of the paper arXiv:2203.16160v1, 
   #         we consider throughout a circuit 
   #         which consists of three blocs) 
   #
   #   iii)  initialization : random starting values  
   #         for all variables in the circuit 
   #         taking into account the structure of the circuit 
   #
   #   iv)   definition of euler steps : update 
   #         of single neurons which receive input 
   #         either inhibitory or excitatory 
   #  
   #     v)  simulation of the circuit as a function of time 
   #         via euler schemes of step size  delta  : 
   #         a matrix  feld  is used to store all variables as a function of time 
   #         graphics illustrate the spiking patterns around the circuit 
   # 
   #    vi)  graphical representations of spiking patterns 
   #         along the circuit 
   #         of interacting stochastic Hodgkin Huxley neurons 




#########################################
#
#    II i) construction of transfer functions according to 5.4 
#          of the paper arXiv:2203.16160v1 submitted to MNA 
#
#########################################
   
   
   
   # general assumption : above, we did chose  tau  and  sigma   
   # such that the stochastic Hodgkin Huxley neuron 
   # is quiet under signal  a_1  
   # and regularly spiking under signal  a_2 : 
   # as explained in example 5.2 of the paper arXiv:2203.16160v1 submitted to MNA 
   
a1 <- 4 ; # neuron is eventually quiet  
a2 <- 10 ; # neuron is eventually regularly spiking   

   # these a_i  were termed  N0_input_faktor  above 
   # in the paper arXiv:2203.16160 submitted to MNA, 
   # we write  theta_1 und theta_2  instead of  a_1  and a_2  
 
 
 
   # take the  u_inf  and   u_sup   as determined in  I vii)  above 
   # for a neuron under  N0_input_faktor = a_2 = 10  
u_inf ; 
u_sup ; 
   
   # from  a_1 , a_2  and  u_inf , u_sup  define now the transfer functions 
   # as in the example following 5.4 of the paper arXiv:2203.16160 submitted to MNA : 
   
transfer_exc <- function( x ){
   m <- ( u_inf + 1 )/2 ; 
   s <- ( u_inf - 1 )/6  ; 
   a1 + (a2-a1) * pnorm( x, m, s ) ; 
   } ; # the exciting transfer function 

transfer_inh <- function( x ){
   m <- ( u_inf + 1 )/2 ; 
   s <- ( u_inf - 1 )/6  ; 
   a1 + (a2-a1) * ( 1 - pnorm( x, m, s ) ) ; 
   } ; # the inhibiting transfer function 
  
  
  
   # graphical representation of the transfer functions : 
leg <- "the neuron is quiet under input  a1 = %%%% and regularly spiking under input  a2 = &&&& ; under input  a2, the  output  u  will eventually fluctuate between   u_inf = §§§§ und u_sup = !!!!  (benchmarks)" ; 
leg <- gsub( "§§§§", round(u_inf,2), leg ) ; 
leg <- gsub( "!!!!", round(u_sup,2), leg ) ; 
leg <- gsub( "%%%%", round(a1,1), leg ) ; 
leg <- gsub( "&&&&", round(a2,1), leg ) ; leg ; 
#
leg1 <- " excitation : transfer function 'Phi_exc' " ;  
leg2 <- " inhibition : transfer function 'Phi_inh' " ; 
leg3 <- "the (0,1)-valued control function" ; 
#
par(mfrow=c(2,1) ) ; 
#
xhilf <- seq( 0 , round(u_sup+1,2) , 0.01 ) ; 
yhilf <- transfer_exc( xhilf  ) ; 
plot( range(xhilf), range(yhilf), type="n", main=leg1, xlab=leg, ylab="" ) ;  
lines( xhilf, yhilf, lwd=2, col=2 ) ; 
abline( v = u_inf, lty=3, col=4 ) ; 
abline( v = u_sup, lty=3, col=4 ) ; 
abline( v = 1, lty=3 ) ;  
#
yhilf <- transfer_inh( xhilf ) ;   
plot( range(xhilf), range(yhilf), type="n", main=leg2, xlab=leg, ylab="" ) ;  
lines( xhilf, yhilf, lwd=2, col=4 ) ; 
abline( v = u_inf, lty=3, col=4 ) ; 
abline( v = u_sup, lty=3, col=4 ) ; 
abline( v = 1, lty=3 ) ;  
#
par(mfrow=c(1,1) ) ; 




#############################################
#
#   II ii) fixing the size of the circuit, and other preparations : 
#
#          as in to 5.5 of the paper arXiv:2203.16160v1 submitted to MNA 
#          except that we spezialize throughout to a circuit with  3  blocs 
#
#          we prepare vectors which will be used in the sequel 
#          to store all variables for all neurons as a function of time  
#
#############################################


     # it is important that all previous parts of the programm 
     # have been performed before passing to  II ii) 

   
     #  we shall define a circuit of  M = 3 L  neurons
     #  counting modulo  M   from  0  to  M-1 ; the transfers     
     #
     #  neuron 0 --> neuron 1 , neuron L --> L + 1 , ... , neuron 2 L --> 2 L + 1  
     #
     #  will be inhibiting, and all other transfers exciting 
     
     
     #  we shall define the matrix  feld  to store 
     #  all variables for all neurons as a function of time  
     #
     #  in the matrix  feld , we reserve the rows 
     #
     #  (i-1)*10 + (1:10) 
     #
     #  for neuron  i 



   #  verification of the setting  : 
a1 ;  # we write  theta_1  instead of  a_1  in the paper arXiv:2203.16160 submitted to MNA 
a2 ;  # we write  theta_2  instead of  a_2  in the paper arXiv:2203.16160 submitted to MNA 
delta ; # the time steps of the euler schemes 
u_inf ;   # lower benchmark for the output of the regularly spiking neuron in the long run 
u_sup ;   # upper benchmark for the output of the regularly spiking neuron in the long run 
lambda ; # we write  c_1  instead of  lambda  in the paper arXiv:2203.16160 submitted to MNA 
   # one has to check if choice of  lambda  is compatible with  5.3 in the paper 
 
 
 
    #  size of the circuit : 
L <- 4 ;  
# L <- 3 ; 
# L <- 7 ; 
L ; # the bloc length 
M <- 3 * L ; # we always work with 3 blocs 
M ;  # total length of the circuit  

   # in the paper,we write  N  for the total length of the circuit 
 
 
 
    # the time horizon and the time grid : 
# T_lim <- 600 ;  # better to start carefully 
T_lim <- 1200 ; 
# T_lim <- 1800 ;  # this (or more) perhaps at the end 
tt <- seq( 0, T_lim, delta ) ; 
length(tt) ; 
   
  
  
  
  #  introduce a matrix  feld  
  #  to store the trajectories of all variables of all neurons   
  
feld <- matrix( 0, nrow=10*M, ncol=length(tt) ) ; 
dim(feld) ; 
feld[ ,1 ] ; 

  # where we attribute  10  rows to every neuron 
  # such that the rows 
  #
  #  (i-1)*10 + (1:10) 
  #
  # are reserved for neuron  i  
   
  # for neuron 1 , 
  # lines 1 , .. , 10 are attributed as follows : 
  # feld[ 1, ] <- N0_input ; 
  # feld[ 2, ] <- N0_phi ; 
  # feld[ 3, ] <- N0_v ; 
  # feld[ 4, ] <- N0_n ; 
  # feld[ 5, ] <- N0_m ;
  # feld[ 6, ] <- N0_h ;
  # feld[ 7, ] <- N0_spike ;
  # feld[ 8, ] <- N0_output ; 
  # and the last two rows will remain unattributed and filled with  0 
  
  # analogous notation modulo 10 for all other neurons : 
  # feld[ (i-1)*10 + 1, ] <- input for neuron  i   
  # ....  
  # feld[ (i-1)*10 + 8, ] <- output generated by neuron  i  
  # with again the last two rows will remain unattributed and filled with  0  




#####################################################
#
#   II iii) initialization : choice of starting values  
# 
#           suitably at random for the Ornstein Uhlenbeck processes, 
#           the biological variables and the output processes 
#
#           in accordance with the structure of the circuit 
#           for the input variables 
#
#####################################################


   # in a first step, we fix at random the starting values 
   # for membrane potential, gating variables and output : 
   
   # random from invariant measure for the Ornstein Uhlenbeck processes 
   # random on (-12,120) for the membrane potential, 
   # random on (0,1) for the gating variables, 
   # random on (1,u_inf) for the output variables (or, alternatively : deterministically = 0 always)  
   
   # in a second step, we initialize the input variables according to section 5 of the paper 
   # as a function of the output of their predecessor 
   # and of its position in the circuit (i.e.: predecessor in the same bloc or not ?) 
   
   
   
   # a first function deals only with the single neuron 
   # and initializes all variables except the input variable : 
   
initialisiere_einzelneuron <- function( ){ 
   startwert <- rep(0,10) ; # preliminary attribution 
   startwert[1] <- runif( 1, a1, a2 ) ; # a preliminary attribution which 
   # has to be redefined according to the output of the predecessor 
   # and the position of the neuron in the circuit 
   startwert[2] <- rnorm( 1, 0, gamma/sqrt(2) )  ; # sampling from the invariant law  
   # of the Ornstein Uhlenbeck process  phi 
   startwert[3] <- runif( 1, -12, 120 )  ; # variable  v
   startwert[4] <- runif( 1, 0, 1 )  ; # variable  n
   startwert[5] <- runif( 1, 0, 1 )  ; # variable  m
   startwert[6] <- runif( 1, 0, 1 )  ; # variable  h 
   startwert[7] <- 0 ; # impossible to have a spike at time  0 
   startwert[8] <- runif( 1, 1, u_inf )  ; # variable u  output 
   # alternative scenario : replace the last line by 
   # startwert[8] <- 0 ;  # output process starts in  0  
   startwert ; 
   } ; # the function initializes all variables of a single neuron except the input 
   # (which has to depend on the predecessor and on the structure of the circuit) 
   # test : 
   # neuron_alt <- rep(0,10) ; 
   # neuron_neu <- initialisiere_einzelneuron( ) ; 
   # neuron_neu ; 



   # a second function initializes all neurons in the circuit 
   # except for the input variables (which have to be defined later 
   # as a function of the output of the predecessor 
   # and of its position in the circuit) : 
   
initialisiere <- function( M ){ 
   hilf <- rep( 0, M*10 ) ; 
   for( i in 1:M ) hilf[ (i-1)*10 + (1:10) ] <- initialisiere_einzelneuron( ) ; 
   hilf ; 
   }; # the function initializes all neurons in the circuit 
   # except the input  
   # test : 
   # feld[ ,1] ; # first row corresponds to time  0 
   # feld[ ,1 ] <- initialisiere( M ) ; 
   # feld[ ,1 ] ; 
   # for( i in 0:(M-1) ) cat( "\n", "neuron ", i+1, ":", "\t", feld[ i*10 + (1:10) , 1 ], "\n" ) ;  
   
   
   
   # so far, the starting values of the input variables 
   # are not yet functions of the output of the preceding neuron 
   # and of its position in the circuit   
   # as required in section 5 of the paper arXiv:2203.16160 submitted to MNA 


 
   # now redefine the input variables correctly 
   # depending on the output of the preceding neuron 
   # and its position in the circuit  (i.e., within the same bloc or not ?) 
   
adaptiere <- function( circuit_startwert ){ 
   hilf <- circuit_startwert ; 
   # the inhibitory transfers :
   hilf[ 0*10 + 1 ] <- transfer_inh( hilf[ (M-1)*10 + 8 ] ) ; # input for neuron 1 
   hilf[ L*10 + 1 ] <- transfer_inh( hilf[ (L-1)*10 + 8 ] ) ; # input for neuron L+1 
   hilf[ 2*L*10 + 1 ] <- transfer_inh( hilf[ (2*L-1)*10 + 8 ] ) ; # input for neuron 2L+1 
   # the excitatory transfers : 
   for( i in 1:(L-1) ) hilf[ i*10 + 1 ] <- transfer_exc( hilf[ (i-1)*10 + 8 ]  ) ;  # input for neuron i+1 
   for( i in (L+1):(2*L-1) ) hilf[ i*10 + 1 ] <- transfer_exc( hilf[ (i-1)*10 + 8 ]  ) ;  # input for neuron i+1 
   for( i in (2*L+1):(M-1) ) hilf[ i*10 + 1 ] <- transfer_exc( hilf[ (i-1)*10 + 8 ]  ) ;  # input for neuron i+1 
   hilf ; 
   } ; # function modifies the preliminary definition of input for all neurons in the circuit
   # input is now a function of the output of the preceding neuron   
   # neurons  1 , L+1 , 2L+1  (those in first position of their bloc) 
   # are inhibited by their predecessor 
   # all other neurons (those having their predecessor inside the same bloc) 
   # are excited by their predecessor 





   # test for all three functions  
   # initialisiere_einzelneuron ,  initialisiere ,  adaptiere  
   # which have been defined above   
   
feld[ ,1 ] <- initialisiere( M ) ; 
for( i in 0:(M-1) ) cat( "\n", "starting value (preliminary) for neuron ", i+1, ":", "\t", feld[ i*10 + (1:10), 1 ], "\n" ) ;  
   # 
circuit_startwert <- feld[ ,1 ] ; 
feld[ ,1 ] <- adaptiere( circuit_startwert  ) ; 
for( i in 0:(M-1) ) cat( "\n", "starting value (adapted) for neuron ", i+1, ":", "\t", feld[ i*10 + (1:10), 1 ], "\n" ) ;  



   # now the first column  feld[,1]   of the matrix  feld   
   # contains random starting values for 
   # the biological variables
   # the driving OU processes, 
   # the output variables 
   # whereas input variables are defined as functions of the output of the preceding neuron 
   # spike indicators are always  0  since there can be no spike at time  0  
   # this is our set of starting values for the euler scheme 
   



 
####################################
#
#  II iv) : prepare functions in order to perform 
#           one euler step for one particular neuron 
#           (fixed time step  delta ) : 
#           the neuron receives either inhibitory or excitatory input;    
#           for this particular neuron, all variables will be updated 
#
#####################################
 
 
    # we introduce two functions : 
    # the first function updates 
    # neurons in the circuit which receive inhibitory input 
    # the second function updates 
    # neurons in the circuit which receive excitatory input 
    
 
   # inhibition : 
   
   # euler stpe update (the time step is  delta  defined above, and now fixed): 
   # new state for neuron i  
   # when transfer from  i-1 (modulo M) to  i  is inhibitory : 
   
eulerschritt_inh <- function( neuron_alt, u_alt_vorgaenger ){ 
   # 'alt=old' refers to the system before performing the euler step 
   neuron_neu <- rep(0,10) ; 
   #  'neu=new' refers to the system after the euler step is performed 
   #
   input <- transfer_inh( u_alt_vorgaenger ) ;  
   # the input is a function of the output of the preceding neuron 
   phi_alt <- neuron_alt[2] ;  # ornstein uhlenbeck process phi   
   v_alt <- neuron_alt[3] ;   # the biological variables 
   n_alt <- neuron_alt[4] ; 
   m_alt <- neuron_alt[5] ; 
   h_alt <- neuron_alt[6] ; 
   # the 'old' spiking pattern variable is not needed 
   output_alt <- neuron_alt[8] ;  # the output  u   
   # the lines  9 und 10  for every neuron (filled with  0 ) are not attributed 
   # 
   # now the updates 
   dphi <- d_phi( phi_alt, delta, gamma, tau ) ; # ornstein uhlenbeck increment 
   dif <- d_IF( v_alt, n_alt, m_alt, h_alt, delta ) ; 
   v_neu <- v_alt + ( input*delta + dphi ) - dif ;  # membrane potential 
   n_neu <- n_alt + d_n( v_alt, n_alt, delta ) ; 
   m_neu <- m_alt + d_m( v_alt, m_alt, delta ) ; 
   h_neu <- h_alt + d_h( v_alt, h_alt, delta ) ;  # gating variables 
   dN <- ifelse( m_alt<h_alt && m_neu>h_neu, 1, 0 ) ; # spiking pattern 
   output_neu <-  output_alt - lambda*output_alt*delta + c*dN  ;  # output 
   #
   neuron_neu[1] <- input ; 
   neuron_neu[2] <- phi_alt + dphi ; 
   neuron_neu[3] <- v_neu ; 
   neuron_neu[4] <- n_neu ; 
   neuron_neu[5] <- m_neu ; 
   neuron_neu[6] <- h_neu ; 
   neuron_neu[7] <- dN ; 
   neuron_neu[8] <- output_neu  ; 
   neuron_neu[9] <- 0 ; 
   neuron_neu[10] <- 0 ; 
   #
   neuron_neu ; 
   } ; # function calculates one euler step : 
   # update over a time interval of length  delta  
   # for a neuron which receives inhibitory input 




   # excitation : 
   
   # euler step update (the time step is  delta  defined above, and now fixed): 
   # new state for neuron i    
   # when transfer from  i-1 (modulo M) to  i  is excitatory  
   
eulerschritt_exc <- function( neuron_alt, u_alt_vorgaenger ){ 
   # 'alt=old' refers to the system before performing the euler step 
   neuron_neu <- rep(0,10) ; 
   # 'neu=new' refers to the system after the euler step is performed 
   #
   input <- transfer_exc( u_alt_vorgaenger ) ;  
   # the input is a function of the output of the preceding neuron 
   phi_alt <- neuron_alt[2] ;  # ornstein uhlenbeck process phi  
   v_alt <- neuron_alt[3] ;   # the biological variables 
   n_alt <- neuron_alt[4] ; 
   m_alt <- neuron_alt[5] ; 
   h_alt <- neuron_alt[6] ; 
   # the 'old' spiking pattern variable is not needed 
   output_alt <- neuron_alt[8] ;  # output variable  u  
   # the lines  9 und 10  for every neuron (filled with  0 ) are not attributed 
   # 
   # now the updates 
   dphi <- d_phi( phi_alt, delta, gamma, tau ) ; # ornstein uhlenbeck increment 
   dif <- d_IF( v_alt, n_alt, m_alt, h_alt, delta ) ; 
   v_neu <- v_alt + ( input*delta + dphi ) - dif ; # membrane potential 
   n_neu <- n_alt + d_n( v_alt, n_alt, delta ) ; 
   m_neu <- m_alt + d_m( v_alt, m_alt, delta ) ; 
   h_neu <- h_alt + d_h( v_alt, h_alt, delta ) ; # gating variables 
   dN <- ifelse( m_alt<h_alt && m_neu>h_neu, 1, 0 ) ; # spiking pattern 
   output_neu <-  output_alt - lambda*output_alt*delta + c*dN    ; # output  u 
   #
   neuron_neu[1] <- input ; 
   neuron_neu[2] <- phi_alt + dphi ; 
   neuron_neu[3] <- v_neu ; 
   neuron_neu[4] <- n_neu ; 
   neuron_neu[5] <- m_neu ; 
   neuron_neu[6] <- h_neu ; 
   neuron_neu[7] <- dN ; 
   neuron_neu[8] <- output_neu  ; 
   neuron_neu[9] <- 0 ; 
   neuron_neu[10] <- 0 ; 
   #
   neuron_neu ; 
   } ; # function calculates one euler step : 
   # update over a time interval of length  delta  
   # for a neuron which receives excitatory input  

   
 


   # test : 
   # from state of neuron  i  ( rows (i-1)*10 + 1:10 )  at time  0*delta  (column 1) 
   # calculate state for neuron i  ( rows (i-1)*10 + 1:10 ) at time  1*delta  (column 2) : 
   # with identification  neuron 0  =  neuron M  
i <- 1 ; 
# i <- L+1 ; # transfer to neuron  L+1  is inhibitory 
# i <- L ; # transfer to neuron  L  is inhibitory 
prec_nr <- ifelse( i==1 , M, i-1 ) ;  # the predecessor of neuron i 
i ; 
prec_nr ; 

neuron_alt <- feld[ (i-1)*10 + (1:10), 1 ] ; # state of neuron  i  at time  0*delta 
  # column  j=1  corresponds to time  (j-1)*delta = 0*delta 
neuron_alt ; 

neuron_vorgaenger <- feld[ (prec_nr-1)*10 + (1:10), 1 ] ; # state of the predecessor neuron at time  0*delta   
   # column  j=1  corresponds to time  (j-1)*delta = 0*delta 
neuron_vorgaenger ; 
  
u_alt_vorgaenger <- feld[ (prec_nr-1)*10 + 8, 1 ] ;  # output written in position  8  in all neurons 
u_alt_vorgaenger ; 
if( i%%L == 1 )  cat( "\n", "transfer is inhibitory", "\n" ) else cat( "\n", "transfer is excitatory", "\n" ); 

  
neuron_neu  <- eulerschritt_exc( neuron_alt, u_alt_vorgaenger ) ; 
if( i%%L == 1 ) neuron_neu  <- eulerschritt_inh( neuron_alt, u_alt_vorgaenger ) ; 
feld[ (i-1)*10 + (1:10), 2 ] <- neuron_neu ;  
neuron_neu ; 
  # the new state of neuron i at time  1*delta 
  # column  j=2  corresponds to time  (2-1)*delta  
  
  # compare old and new states  
cat( "\n", "neuron ", i, " before euler step :", "\t", neuron_alt, "\n" ) ;  
cat( "\n", "output of neuron ", prec_nr , " which is the predecessor of neuron ", i ," : ", feld[ (prec_nr-1)*10 + 8 , 1 ], "\n" ) ; 
cat( "\n", "neuron ", i, "  after euler step :", "\t", neuron_neu, "\n" ) ;  
   
   
    

  
   # now we are ready to simulate the circuit as a whole 
   
   


#########################################
#
#    II v) simulation of a circuit of  M = 3*L  neurons 
#
#    the simulation run 
#
#    we obtain the data 
#    needed to draw pictures of spiking patterns comparable to  figures 7 or 8 
#    in the paper arXiv:2203.16160v1 submitted to MNA 
#
#########################################


 
   # time horizon for the circuit (better to start carefully) : 
# T_lim <- 600 ;   
T_lim <- 1200 ; 
# T_lim <- 1800 ; 
# T_lim<- 2400 ; 
tt <- seq( 0, T_lim, delta ) ; 
length(tt) ; 
  # the time grid 
  
  
  # the matrix  feld  
  # where we shall store all variables of all neurons as a function of time 
feld <- matrix( 0, nrow=10*M, ncol=length(tt) ) ; 
dim(feld) ; 
feld[ ,1 ] ; 


   # initialization 
circuit_startwert <- initialisiere( M ); 
   # these are the random starting values 
   # all entries  input  are still preliminary and have to be adapted 
feld[,1] <- adaptiere( circuit_startwert ) ; 
feld[,1] ; 
for( i in 1:M ) cat( "\n", "neuron ", i, " starts in state " , "\t", feld[ (i-1)*10 + (1:10), 1 ], "\n" ) ; 
   # initialization of the circuit is now finished 
   # and input (inhibitory / excitatory) for a successor neuron is now 
   # function of the output of the predecessor 
   # according to the structure of the circuit 


   # now the starting configuration is stored in the first column of  feld  


   # simulation of the circuit as a whole : 
   # we calculate the matrix  feld  
   # (column  k  corresponds to time  (k-1)*delta , 
   # rows   (i-1)*10 + (1:10)   correspond to neuron  i ) 
   
j <- 2 ; while( j <= length(tt) ){ 
   hilf_alt <- feld[ , j-1 ] ; # state of the circuit at time  (j-1-1)*delta 
   hilf_neu <- feld[ , j ] ; # here we shall note the state of the circuit at time  (j-1)*delta 
   # hilf_neu ; 
   neuron_alt <- rep( 0, M ) ; # preliminary : store the 'old' values of the neuron 
   neuron_neu <- rep( 0, M ) ; # preliminary : store the 'new' values of the neuron 
   #
   # we calculate first the neurons receiving inhibitory transfer 
   # i = 1 
   neuron_alt <- hilf_alt[ 0*10 + (1:10) ] ; # 'old' state of neuron 1 
   u_alt_vorgaenger <- hilf_alt[ (M-1)*10 + 8 ] ; # output of the preceding neuron  M = 0 
   neuron_neu <- eulerschritt_inh( neuron_alt, u_alt_vorgaenger ) ; 
   hilf_neu[ 0*10 + (1:10) ] <- neuron_neu ;  # 'new' state of neuron 1   
   # i = L+1 
   neuron_alt <- hilf_alt[ L*10 + (1:10) ] ; # 'old' state of neuron  i = L+1
   u_alt_vorgaenger <- hilf_alt[ (L-1)*10 + 8 ] ; # output of the preceding neuron  i-1 = L 
   neuron_neu <- eulerschritt_inh( neuron_alt, u_alt_vorgaenger ) ; 
   hilf_neu[ L*10 + (1:10) ] <- neuron_neu ;  # 'new' state of neuron  i = L+1 
   # i = 2L+1 
   neuron_alt <- hilf_alt[ (2*L)*10 + (1:10) ] ; # 'old' state of neuron i = 2L+1 
   u_alt_vorgaenger <- hilf_alt[ (2*L-1)*10 + 8 ] ; # output of the preceding neuron  i-1 = 2L
   neuron_neu <- eulerschritt_inh( neuron_alt, u_alt_vorgaenger ) ; 
   hilf_neu[ (2*L)*10 + (1:10) ] <- neuron_neu ;  # 'new' state of neuron i = 2L+1 
   #
   # now we calculate the neurons receiving excitatory transfer  : 
   for( i in 2:L ){  
   neuron_alt <- hilf_alt[ (i-1)*10 + (1:10) ] ; # 'old' state of neuron  i 
   u_alt_vorgaenger <- hilf_alt[ (i-1-1)*10 + 8 ] ; # output of the preceding neuron  i-1 
   neuron_neu <- eulerschritt_exc( neuron_alt, u_alt_vorgaenger ) ; 
   hilf_neu[ (i-1)*10 + (1:10)  ] <- neuron_neu ;  # 'new' state of neuron i
   } ;   
   for( i in (L+2):(2*L) ){
   neuron_alt <- hilf_alt[ (i-1)*10 + (1:10) ] ; # 'old' state of neuron  i 
   u_alt_vorgaenger <- hilf_alt[ (i-1-1)*10 + 8 ] ; # output of the preceding neuron  i-1 
   neuron_neu <- eulerschritt_exc( neuron_alt, u_alt_vorgaenger ) ; 
   hilf_neu[ (i-1)*10 + (1:10)  ] <- neuron_neu ; # 'new' state of neuron i
   } ;  
   for( i in (2*L+2):M ) {
   neuron_alt <- hilf_alt[ (i-1)*10 + (1:10) ] ; # 'old' state of neuron  i 
   u_alt_vorgaenger <- hilf_alt[ (i-1-1)*10 + 8 ] ; # output of the preceding neuron  i-1 
   neuron_neu <- eulerschritt_exc( neuron_alt, u_alt_vorgaenger ) ; 
   hilf_neu[ (i-1)*10 + (1:10)  ] <- neuron_neu ; # 'new' state of neuron i
   } ;  
   # 
   feld[ , j ] <- hilf_neu ;  # 'new' state of the whole circuit 
   #
   cat( "\n", "time ", round( (j-1)*delta, 2 ), "\n" ) ; 
   j <- j+1 ; 
   } ; # ende while  
   # this loop calculates the circuit up to time   T_lim  
   # and stores all variables of all variables as a function of time 
   # in the matrix  feld  



   # now the circuit of neurons as a whole is calculated 
   # and stored in the matrix  feld 
   # (column  k  corresponds to time  (k-1)*delta , 
   # rows   (i-1)*10 + (1:10)   correspond to neuron  i ) 
   
   
   
   
###########################################
#
#   II vi) graphical representations of spiking patterns 
#          along the circuit of neurons  
#          calculated in  II v)  above  
#
###########################################
   
   
   # first type of graphics : 
   # one figure showing the spiking patterns around the circuit   
   # (cf. figures 7 and 8 in the paper arXiv:2203.16160v1 submitted to MNA)  
   
leg <- "spike trains in circuit of 3 x §§§§ neurons, inhibited (red) or excited (green) by predecessor (OU sigma = %%%% , tau = !!!! )" ; 
leg <- gsub( "§§§§", L, leg ) ;
leg <- gsub( "%%%%", round(sigma,2), leg ) ; 
leg <- gsub( "!!!!", round(tau,2), leg ) ; 
leg ; 
yleg <- "neuron i on level i, counting modulo §§§§ " ; 
yleg <- gsub( "§§§§", M, yleg ) ; 
yleg ; 
xleg <- "time (from 0 to !!!!, step size %%%%)" ; 
xleg <- gsub( "%%%%", delta, xleg ) ; 
xleg <- gsub( "!!!!", T_lim, xleg ) ; 
xleg ; 
plot( range(tt), c(-2,M+2), type="n", ylab=yleg, xlab="time", main=leg ) ; 
for ( j in 0:M ) abline( h=j, lty=3, cex=0.1 ) ; 
   # neuron  i  on level i , with level M = 3*L shown twice 
   # (i.e. on level  0  and also on level  M ) to emphasize circular structure 
tt0 <- tt[ feld[ (M-1)*10 + 7 , ] == 1 ] ; 
points( tt0, rep( M, length(tt0) ) , col=3, pch=19 ) ;  
points( tt0, rep( 0, length(tt0) ), col=3 ) ; 
   # transfer to neurons (red) at the beginning of the three blocs is inhibitory : 
tt0 <- tt[ feld[ 0*10 + 7 , ] == 1 ] ; 
yy0 <- rep( 0+1, length(tt0) ) ; 
points( tt0, yy0, col=2, pch=19 ) ;  # spiking pattern neuron 1 
tt0 <- tt[ feld[ L*10 + 7 , ] == 1 ] ; 
yy0 <- rep( L+1, length(tt0) ) ; 
points( tt0, yy0, col=2, pch=19 ) ; # spiking pattern neuron L+1  
tt0 <- tt[ feld[ 2*L*10 + 7 , ] == 1 ] ; 
yy0 <- rep( 2*L+1, length(tt0) ) ; 
points( tt0, yy0, col=2, pch=19 ) ;  # spiking pattern neuron 2L+1  
   # while transfer to neurons (green) inside every bloc is excitatory : 
for( k in 2:L){
   tt0 <- tt[ feld[ (k-1)*10 + 7 , ] == 1 ] ; 
   yy0 <- rep( k, length(tt0) ) ; 
   points( tt0, yy0, col=3, pch=19 ) ; 
   } ; # spike patterns inside bloc 1
for( k in (L+2):(2*L) ){
   tt0 <- tt[ feld[ (k-1)*10 + 7 , ] == 1 ] ; 
   yy0 <- rep( k, length(tt0) ) ; 
   points( tt0, yy0, col=3, pch=19 ) ; 
   } ;  # spike patterns inside bloc 2
for( k in (2*L+2):M){
   tt0 <- tt[ feld[ (k-1)*10 + 7 , ] == 1 ] ; 
   yy0 <- rep( k, length(tt0) ) ; 
   points( tt0, yy0, col=3, pch=19 ) ; 
   } ;  # spike patterns inside bloc 3
   # now the figure is complete  




  # second type of graphics : 
  # a series of figures showing   
  #  1) output  U  delivered by its predecessor 
  #  2) input  A  received by the neuron (red color if 
  #     neuron is inhibited by predecessor, green if excited)  
  #  3) membrane potential  V  
  # for all variables in the circuit as a function of time 
  
  

  # a first slide shows neuron 0 == M = 3*L  : 
leg3 <- "circuit of  §§§§  neurons arranged in 3 blocs of length  %%%%  :  output  U  delivered by neuron §§§§ - 1 , the predecessor of neuron  0 == §§§§" ; 
leg3 <- gsub( "§§§§", M, leg3 ) ; 
leg3 <- gsub( "%%%%", L, leg3 ) ;  
#  
leg1 <- "which transforms into input  A  received by neuron  0 == §§§§" ; 
leg1 <- gsub( "§§§§", M, leg1 ) ; 
#
leg2 <- "membrane potential  V  of neuron  0 == §§§§ , OU parameters : sigma = %%%% , tau = !!!!" ; 
leg2 <- gsub( "§§§§", M, leg2 ) ; 
leg2 <- gsub( "%%%%", round(sigma,2), leg2 ) ; 
leg2 <- gsub( "!!!!", round(tau,2), leg2 ) ; 
#
par(mfrow=c(3,1) ) ;  
plot( range(tt), c(0,1.1*u_sup), type="n", main=leg3, ylab="U", xlab="time" ) ;
abline( h = u_sup, lty=2 ) ; 
abline( h = u_inf, lty=2 ) ; 
uu <- feld[ ((M-1)-1)*10 +8 , ] ;  
lines( tt, uu, col=4, lwd=2 ) ; 
plot( range(tt), c(0.5*a1,1.2*a2), type="n", main=leg1, ylab="A", xlab="time" ) ; 
abline( h = a1, lty=2 ) ; 
abline( h = a2, lty=2 ) ; 
lines( tt, feld[ (M-1)*10 + 1 , ], col=3, lwd=2 ) ;
plot( range(tt), c(-12,120), type="n", main=leg2, ylab="V", xlab="time" ) ;
lines( tt, feld[ (M-1)*10 + 3 , ] ) ; 
par(mfrow=c(1,1) ) ;  
#
  # the next slides show successively the neurons along the circuit : 
i <- 1 ; while( i <= M ){ 
   # 
   col_i <- 3 ;  # color if transfer to neuron  i  is excitatory  
   if( (i%%L) == 1 ) col_i <- 2 ;  # color if transfer is inhibitory   
   #
   leg3 <- "circuit of  %%%%  neurons arranged in 3 blocs of length  !!!!  :  output  U  delivered by neuron  §§§§ - 1 , the predecessor of neuron  §§§§" ; 
   leg3 <- gsub( "§§§§", i, leg3 ) ; 
   leg3 <- gsub( "%%%%", M, leg3 ) ; 
   leg3 <- gsub( "!!!!", L, leg3 ) ; 
   # 
   leg1 <- "transfer from  §§§§ - 1  to  §§§§ is exciting : input  A  received by neuron §§§§" ; 
   if( (i%%L) == 1 ) leg1 <- "transfer from  §§§§ - 1  to  §§§§ is inhibiting : input  A  received by neuron §§§§" ; 
   leg1 <- gsub( "§§§§", i, leg1 ) ;  
   # 
   leg2 <- "membrane potential  V  of neuron §§§§ , OU parameters : sigma = %%%% , tau = !!!!" ; 
   leg2 <- gsub( "§§§§", i, leg2 ) ; 
   leg2 <- gsub( "%%%%", round(sigma,2), leg2 ) ; 
   leg2 <- gsub( "!!!!", round(tau,2), leg2 ) ; 
   # 
   par( mfrow=c(3,1) ) ; 
   plot( range(tt), c(0,1.1*u_sup), type="n", main=leg3, ylab="U", xlab="time" ) ;
   abline( h = u_sup, lty=2 ) ; 
   abline( h = u_inf, lty=2 ) ; 
   if( i == 1 ) uu <- feld[ (M-1)*10 +8 , ] ; 
   if( i > 1 ) uu <- feld[ ((i-1)-1)*10 +8 , ]  ;  
   lines( tt, uu, col=4, lwd=2 ) ; 
   plot( range(tt), c(0.5*a1,1.2*a2), type="n", main=leg1, ylab="A", xlab="time" ) ; 
   abline( h = a1, lty=2 ) ; 
   abline( h = a2, lty=2 ) ; 
   lines( tt, feld[ (i-1)*10 + 1 , ], col=col_i, lwd=2 ) ;
   plot( range(tt), c(-12,120), type="n", main=leg2, ylab="V", xlab="time" ) ;
   lines( tt, feld[ (i-1)*10 + 3 , ] ) ; 
   par( mfrow=c(1,1) ) ;  
   #
   i <- i+1 ; 
   } ;  # end of loop 




   # a third type of graphics : 
   # membrane potential of all neurons in a bloc 
   # around the circuit 
   # beginning and ending with neuron 0 == M  

leg0 <- "circuit of §§§§ neurons arranged in 3 blocs of length  %%%% , spike train of neuron 0 = neuron §§§§" ; 
leg0 <- gsub( "§§§§", M, leg0 ) ; 
leg0 <- gsub( "%%%%", L, leg0 ) ; 
plot( tt, feld[ (M-1)*10 +3, ], type="l", main=leg0, ylab="V", col=2 ) ; 
par( mfrow=c(L,1) ) ; 
for( m in 1:M ){
   leg <- "circuit of §§§§ neurons arranged in 3 blocs of length  %%%% , spike train of neuron !!!!" ; 
   leg <- gsub( "!!!!", m, leg ) ; 
   leg <- gsub( "§§§§", M, leg ) ; 
   leg <- gsub( "%%%%", L, leg ) ; 
   plot( tt, feld[ (m-1)*10 +3, ], type="l", main=leg, ylab="V", col=2 ) ; 
   } ; 
par(mfrow=c(1,1) ) ;  
plot( tt, feld[ (M-1)*10 +3, ], type="l", main=leg0, ylab="V", col=2 ) ; 





#############################################
#
#   end of program 
#
#   last update  26.08.22
#
#############################################





