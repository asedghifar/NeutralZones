#Author: Alisa Sedghifar
#This file has been submitted as a supplemental file to the manuscript "The Spatial Mixing of Genomes in Secondary Contact Zones", 2015. The scripts here provide a very basic framework for fitting a spatial diffusion model of admixture to weighted LD. CAUTION: These scripts come with no guarantees, and there may be bugs.

#You will need to load the following:
	
	#geo_locations: A vector of the geographic positions of populations
	
	#mean_LD: A list of vectors, each corresponding to a population. Each entry of the vector should be the weighted LD corresponding to the genetic distance in genetic_positions (see below). We used the output of ALDER to obtain these values.
	
	#variance_LD: A list of vectors, each corresponding to a population, as above. Each entry of the vector should be the variance of the weighted LD corresponding to the genetic distances in genetic_positions. We took the variance of the jackknife output of ALDER to obtain these values. Note that this should be an estimate of the variance for the jackknife, and not the sample variance. 
	
	#genetic_positions: A list of vectors, each corresponding to a population, as above. Each entry of the vector should be the genetic positions for which weighted LD has been computed, in centiMorgans. ALDER provides this information in its output.
	


#The following function needs to be integrated over T on (0,1) to get expected LD. The variables are (T: time of recombination), (SIGMA: variance parameter of dispersal), (TAU: time since secondary contact), (R: genetic distance of two markers), (L: geographic position of population)

LD_int= function(T,SIGMA,TAU,R,L){
			exp(-R*T*TAU)*exp((T-1)*(L)^2/(SIGMA^2*TAU*(1-T^2)))/(2*pi*sqrt(1-T^2))
		}


#To compute the sum of squares for all populations:
LD_all = function(V){
	sumsq = sum(
				sapply(1:length(genetic_positions),
					function(i){
						sum((mean_LD[[i]]-exp(V[3])*sapply(genetic_positions[[i]]/100,function(R){integrate(LD_int,0,1,SIGMA=exp(V[2]),TAU=exp(V[1]),R=R,L=geo_locations[i]-V[4])$value}))^2/variance_LD[[i]])	
					}
				)
			)
	return(sumsq)
}


#V is the vector c(log(TAU),log(SIGMA),log(F),L). F is a scaling parameter that reflects the differentiation between the parental populations.
#One way to find V is to use the R nlm function. The starting values should be informed by the knowledge of your system. e.g :

mlV =nlm(LD_all,V<-c(log(200),log(0.1),log(0.001),exp(120)))

#If nlm does not give consistent results, an array could be iterated over:

rates = seq() #The values of SIGMA to search over
times = seq() #The values of TAU to search over
fparam = seq() #The values of F to search over

LD_array = array(dim=c(length(times),length(rates),length(fparam)))

for(k in 1:length(fparam)){
	for(j in 1:length(rates)){
		for(i in 1:length(times)){
	LD_array[i,j,k] = JKV(log(c(times[j],rates[i],fparam[k],exp(124.9))))	
		}	
	}
}

#In the above example, we have specificied a value for L based on the fit to the frequency cline. This could also be made into a fourth dimension if so desired.

#To get the profile likelihood surface:
LD_proflik=apply(LD_array,1,function(x){apply(x,1,min)})

#Remember to take the additive inverse of this to get a surface akin to the likelihood surface. 
