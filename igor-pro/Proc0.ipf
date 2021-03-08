#pragma rtGlobals=1		// Use modern global access method.

function KuramotoController()
variable dimension
variable mode, g, epsilon
variable timemetric, clustersize, N, dt, tinitial
variable  Z, drivingfrequency, tmax, spreadinomega
variable timemetric2, epsiloncounter

dimension=1 //Choose what dimension we are working in (options=1 or 2)
if (dimension ==1)
epsiloncounter=0
endif

//all parameters
g=1 //The initial coupling constant
N=100 //Either the number of oscillators (1D) or the height/width of the square array
//epsilon=1
//timemetric=7 //Zannette's "T" used only for initial conditions if timemetric2 != it. 
//timemetric2=timemetric 
dt=.01 //the "amount of time" for each discrete time step. 
tinitial=1000 //The amount of time that the model runs initial conditions

mode=0//Sets initial conditions: 0->random phases, 1->same phase cluster, 2->plane wave cluster, 3->circles cluster
clustersize=N //cluster=N for whole array to be driven
drivingfrequency = 1 //frequency that the oscillaotrs in the cluster will be forced to move at
spreadinomega = 0 //the standard deviation in the natural frequencies
tmax=20000 //the number of total time steps we are simulating (counts the initial conditions)

init(N, tmax, dimension) //function to initialize the memory to store all values. 

for (timemetric=0; timemetric<11;timemetric=timemetric+1)
	for (epsilon=0; epsilon<.1; epsilon=epsilon+.01)
timemetric2=timemetric


//Initial Condition Functions
connectivityanddelays(g, timemetric, N, dt, spreadinomega, dimension) //calls function to assign connectivity, weights and timedelays
initialconditions0(N, dimension)
if (mode == 0)
clustersize=0 //for m=0 we assume they will run at their natural frequencies
endif

if (mode == 1)
initialconditions1(N, clustersize, dimension)
endif

if (mode==2)
initialconditions2(N, clustersize, dimension)
endif

if (mode==3)
initialconditions3(N, clustersize, dimension)
endif

//if we want a different timemetric than the one we used to set up the intial conditions we change timemetric2
if (timemetric2 != timemetric)
wrongwave(N, dt, timemetric2, dimension) //function to assign new tau's. 
endif

//Now we drive the cluster for the intial time then start the model
driver(dt, drivingfrequency, clustersize, tinitial, N, dimension)
model(dt, tinitial, tmax, N, epsilon, dimension)

//This next part will do some analysis if we are moving through epsilon and timemetric in 1D
if (dimension ==1)
OPs(N, tmax)
OPcompare(tmax)
frequencies(N, tmax)
endif

FinalFreqs(N, dt, tmax, dimension)

ExperimentSave(g, N, timemetric, timemetric2, epsilon, mode, dimension) //Last we save each experiment

//if (dimension==1)
//epsiloncounter=epsiloncounter+1
//endif

endfor
endfor

end

//Creation of Wave Space
function init(N, tmax, dimension)
variable N, tmax, dimension

if (dimension==2)
Make/o/n=(N, N, tmax) phi  //must be big enough to hold largest delay and total run time
Make/o/n=(N, N) phicurrent
Make/o/n=(N^2) phicurrent1D
Make/o/n=(N, N) omega
Make/o/n=(N,N) FinalFrequencies
Make/o/n=(N, N, N, N) connectionmatrix
Make/o/n=(N, N, N, N) alpha
Make/o/n=(N, N, N, N) Kcoupling
Make/o/n=(N, N, N, N) tau
elseif (dimension==1) 
Make/o/n=(N, tmax) phi  //must be big enough to hold largest delay and total run time
Make/o/n=(N) phicurrent
Make/o/n=(N) omega
Make/o/n=(N, N) alpha
Make/o/n=(N,N) connectionmatrix
Make/o/n=(N,N) Kcoupling
Make/o/n=(N,N) tau
Make/o/n=(N) FinalFrequencies
Make/o/n=(4, tmax) OrderParameter, OrderParameter2Cluster, OrderParameter2ClusterConnected
Make/o/n=(N*tmax) tt, freqs
endif

end


function connectivityanddelays(g, timemetric, N, dt, spreadinomega, dimension)
variable g, timemetric, N, dt, spreadinomega, dimension
variable i, j, k, l, Z
wave tau, connectionmatrix, kcoupling, alpha, omega
//This handles the creation of the tau, K, alpha, omega and connection matrices
Z=timemetric/dt
Kcoupling=g //constant coupling
if (dimension ==2 )
for(i=0;i<N;i=i+1)
	for(j=0;j<N;j=j+1)
		for(k=0;k<N;k=k+1)
			for(l=0;l<N;l=l+1)
			tau[i][j][k][l]=  trunc(Z/N*sqrt( (min(abs(i-k), N-abs(i-k)))^2 +(min(abs(j-l), N-abs(j-l)))^2 )) 
			endfor
		endfor
	endfor
endfor
elseif (dimension ==1)
for(i=0;i<=N;i=i+1)
		for(j=0;j<=N;j=j+1)
		tau[i][j]=trunc(Z/N*min(abs(i-j),N-abs(i-j)))
		endfor
endfor
endif
alpha=1
connectionmatrix=1
omega=gnoise(spreadinomega) + 1
end

function initialconditions0(N, dimension)
variable N, dimension
variable i, j, k, l  //indexes of oscillators
wave phi
//this assigns phase randomly 
if (dimension==2)
for(i=0;i<N;i=i+1)
	for(j=0;j<N;j=j+1)
		phi[i][j][0]=enoise(pi)
	endfor
endfor
elseif (dimension==1)
for(i=0;i<N;i=i+1)
	phi[i][0]=enoise(pi)
endfor
endif
end

function initialconditions1(N,clustersize, dimension)
variable N, clustersize, dimension
variable i, j, k, l //indexes of oscillators
wave phi
//this sets all oscillators in the cluster at the same phase
if (dimension==2)
for(i=0;i<clustersize;i=i+1)
	for(j=0;j<clustersize;j=j+1)
		phi[i][j][0]=1
	endfor
endfor
elseif (dimension ==1)
for(i=0;i<clustersize;i=i+1)
		phi[i][0]=1
endfor
endif
end

function initialconditions2(N, clustersize, dimension)
variable N, clustersize, dimension
variable i, j, k, l  //indexes of oscillators
variable currentdelay
wave phi, tau
//this creates plane waves according to the timedelay
if (dimension ==2)
phi[0][0][0]=1
for (j=0;j<clustersize;j=j+1)
	currentdelay=tau[0][0][0][j]
	for(i=0;i<clustersize;i=i+1)
	phi[i][j][0]=phi[0][0][0]+currentdelay
	endfor
endfor
elseif (dimension==1)
phi[0][0]=1
for (j=0;j<clustersize;j=j+1)
	currentdelay=tau[0][j]
	for(i=0;i<clustersize;i=i+1)
	phi[i][0]=phi[0][0]+currentdelay
	endfor
endfor
endif
end


function initialconditions3(N, clustersize, dimension)
variable N, clustersize, dimension
variable t, currentdelay, timewithstep
variable i, j, k, l //indexes of oscillators
wave phi, tau
//this sets every oscillators phase according to the reference
if (dimension==2)
phi[0][0][0]=1
for (j=0;j<clustersize;j=j+1)
	for(i=0;i<clustersize;i=i+1)
	currentdelay=tau[0][0][i][j]
	phi[i][j][0]=phi[0][0][0]+currentdelay
	endfor
endfor
elseif(dimension==1)
phi[0][0]=1
for (j=0;j<clustersize;j=j+1)
	for(i=0;i<clustersize;i=i+1)
	currentdelay=tau[0][j]
	phi[i][0]=phi[0][0]+currentdelay
	endfor
endfor
endif
end

function wrongwave(N, dt, timemetric2, dimension)
variable N, dt, timemetric2, dimension
variable i, j, k, l
variable Z
wave tau
//This changes the taus that will be used in the experiment from the taus that were used in the set up
Z=timemetric2/dt
if(dimension==2)
for(i=0;i<N;i=i+1)
	for(j=0;j<N;j=j+1)
		for(k=0;k<N;k=k+1)
			for(l=0;l<N;l=l+1)
			tau[i][j][k][l]=  trunc(Z/N*sqrt( (min(abs(i-k), N-abs(i-k)))^2 +(min(abs(j-l), N-abs(j-l)))^2 )) 
			endfor
		endfor
	endfor
endfor
elseif (dimension==1)
for(i=0;i<=N;i=i+1)
		for(j=0;j<=N;j=j+1)
		tau[i][j]=trunc(Z/N*min(abs(i-j),N-abs(i-j)))
		endfor
endfor
endif
end

function driver(dt, drivingfrequency, clustersize, tinitial, N, dimension)
variable dt, drivingfrequency, clustersize, tinitial, N, dimension
variable t, timewithstep, i, j
wave phi, omega
//this drives the oscillators that we picked out during intial conditions at a specific frequency
if(dimension==2)
for (t=0;t<=tinitial;t=t+1)
	for(i=0;i<N;i=i+1)
		for(j=0;j<N;j=j+1)
			timewithstep=t+1
    		 		if (i<clustersize && j<clustersize)
    		 		phi[i][j][timewithstep]=phi[i][j][t]+drivingfrequency*dt
    		 		endif
    		 		if (i>=clustersize || j >=clustersize)
    		 		phi[i][j][timewithstep]=phi[i][j][t]+omega[i][j]*dt //0<t<tinitial phase assignments
    		 		endif
		endfor
	endfor
endfor
elseif (dimension==1)
for (t=0;t<=tinitial;t=t+1)
	for(i=0;i<N;i=i+1)
			timewithstep=t+1
    		 		if (i<clustersize)
    		 		phi[i][timewithstep]=phi[i][t]+drivingfrequency*dt
    		 		endif
    		 		if (i>=clustersize)
    		 		phi[i][timewithstep]=phi[i][t]+omega[i][j]*dt //0<t<tinitial phase assignments
    		 		endif
	endfor
endfor
endif
end 

function model(dt, tinitial, tmax, N, epsilon, dimension)
variable  dt, tinitial, tmax,  N, epsilon, dimension
variable t, timewithstep, timewithdelay,  timeadjusted, summation
variable i, j, k, l, Z
wave phi, phicurrent, omega, connectionmatrix, Kcoupling, tau, alpha

print "Started at ", time()
summation=0
if (dimension==2)
//BEGINING OF TIME...finally
for(t=tinitial;t<tmax;t=t+1)
	//Evolution of Phase
	for(i=0;i<N;i=i+1)
		for(j=0;j<N;j=j+1)
		phicurrent[i][j]=phi[i][j][t]
			for(k=0;k<N;k=k+1)
				for(l=0;l<N;l=l+1)
				timewithdelay=t-tau[i][j][k][l]
				summation=summation+Kcoupling[i][j][k][l]*sin(phi[k][l][timewithdelay]-phicurrent[i][j])*connectionmatrix[i][j][k][l]
				endfor
			endfor
		timewithstep=t+1
		phi[i][j][timewithstep]=phicurrent[i][j]+(omega[i][j]+(summation/(N^2)))*dt //N is squared because N is now the “width” of the square and not the actual number of oscillators
		summation=0
		endfor
	endfor
	
	//Evolution of Coupling
	for(i=0;i<N;i=i+1)
		for(j=0;j<N;j=j+1)
			for(k=0;k<N;k=k+1)
				for(l=0;l<N;l=l+1)
				timewithdelay=t-tau[i][j][k][l]
				Kcoupling[i][j][k][l]=Kcoupling[i][j][k][l]+(epsilon*(alpha[i][j][k][l]*cos(phicurrent[i][j]-phi[k][l][timewithdelay])-Kcoupling[i][j][k][l])*dt)
				endfor
			endfor
		endfor
	endfor
	
if (mod(t, 2000) == 0)
print  "completed time, t=", t, " at ", time()
endif

endfor
//END OF TIME
elseif (dimension==1)
//BEGINING OF TIME (1D)
for(t=tinitial;t<=tmax;t=t+1)
			//Evolution of Phase
			for(i=0;i<=N;i=i+1)
			phicurrent[i]=phi[i][t]
				for(j=0;j<=N;j=j+1)
				timewithdelay=t-tau[i][j]
				summation=summation+Kcoupling[i][j]*sin(phi[j][timewithdelay]-phicurrent[i])*connectionmatrix[i][j]
				endfor
			timewithstep=t+1
			phi[i][timewithstep]=phicurrent[i]+(omega[i]+(summation/N))*dt
			summation=0
			endfor
	
			//Evolution of Coupling
			for(i=0;i<=N;i=i+1)  
				for(j=0;j<=N;j=j+1)
				timewithdelay=t-tau[i][j]
				Kcoupling[i][j]=Kcoupling[i][j]+(epsilon*(alpha[i][j]*cos(phicurrent[i]-phi[j][timewithdelay])-Kcoupling[i][j])*dt) 
				endfor
			endfor

	
//if (mod(t, 2000) == 0)
//print  "completed time, t=", t, " at ", time()
//endif
endfor
//END OF TIME
endif
end

function FinalFreqs(N, dt, tmax, dimension)
variable N, dt, tmax, dimension
variable i, j
variable deltat, comparisontime
wave phi, finalfrequencies

//calculates the average frequency each oscillator is moving at by the end
deltat=20
if (dimension==1)
for(i=0;i<N;i=i+1)
	comparisontime=tmax-(deltat+1)
	finalfrequencies[i]=(phi[i][tmax]-phi[i][comparisontime])/deltat
endfor

elseif(dimension==2)
for(i=0;i<N;i=i+1)
	for(j=0;j<N;j=j+1)
		comparisontime=tmax-(deltat+1)
		finalfrequencies[i][j]=(phi[i][j][tmax]-phi[i][j][comparisontime])/deltat
	endfor
endfor
endif

finalfrequencies=finalfrequencies*dt
print "Final Frequency info:" 
WaveStats finalfrequencies
end

function ExperimentSave(g, N, timemetric, timemetric2, epsilon, mode, dimension)
variable g, N, timemetric, timemetric2, epsilon, mode, dimension
string  filename
filename = "N" + num2istr(N)  +"T" + num2istr(timemetric2) + "E" + num2istr(epsilon*100) + "M" + num2istr(mode) + "D" + num2istr(dimension)
SaveExperiment as filename
print "Saved", filename, " at ", time()
end

function playback()
variable t, i, j, q, N
wave phi
wave phicurrent, phicurrent1D, omega
//N=5
for (t=0;t<10000;t=t+2)
	for(i=0;i<21;i=i+1)
		for(j=0;j<21;j=j+1)
		phicurrent[i][j]= mod(phi[i][j][t], 2*pi)
		//phicurrent1D[i][j] = mod(phi[i][j][t], 2*pi); redimension/N=(N^2) phicurrent1D
			if(phicurrent[i][j] > pi)
			phicurrent[i][j]= 2*pi-phicurrent[i][j]
			endif
		endfor
	endfor
//duplicate/O omega, omega1D; redimension/N=(N^2) omega1D
Doupdate
endfor
end

function playback1D()
variable t, i, j,N
wave phi
wave phicurrent, phicurrent1D, omega
//N=5
for (t=0;t<20000;t=t+10)
	for(i=0;i<101;i=i+1)
		phicurrent[i]= mod(phi[i][t], 2*pi)
		//phicurrent1D[i][j] = mod(phi[i][j][t], 2*pi); redimension/N=(N^2) phicurrent1D
			//if(phicurrent[i] > pi)
			//phicurrent[i]= 2*pi-phicurrent[i]
			//endif
		endfor


Doupdate
endfor
end

function makemovie()
variable t, i, j, q
wave phi
wave phicurrent, phicurrent2
NewMovie/L/F=30 
for (t=0;t<10000;t=t+10)
	for(i=0;i<100;i=i+1)
		for(j=0;j<100;j=j+1)
		phicurrent[i][j]= mod(phi[i][j][t], 2*pi)
			if(phicurrent[i][j] > pi)
			phicurrent[i][j]= 2*pi-phicurrent[i][j]
			endif
		endfor
	endfor
Doupdate
addmovieframe 
endfor
closemovie
end




function frequencies(N, tmax)
variable N, tmax
variable timewithstep, t, i
wave freq, tt, phi

for(t=0;t<tmax;t=t+1)
	for(i=0;i<N;i=i+1)
	timewithstep=t+1
	freq[i+t*N]=N*(phi[i][timewithstep]-phi[i][t])
	tt[i+t*N]=t
	endfor
endfor

end



function OPs(N, tmax)
variable N, tmax
variable i,t, sinsum, cossum, sinsum2, cossum2, sinsum3, cossum3, m
variable OPminus, OPplus, OPminus2, OPplus2, OPminus3, OPplus3
wave  OrderParameter, OrderParameter2Cluster, OrderParameter2ClusterConnected, phicurrent, phi

//Order Parameter Calculation 
sinsum=0;cossum=0
sinsum2=0;cossum2=0
sinsum3=0;cossum3=0
OPminus=0;OPplus=0
OPminus2=0;OPplus2=0
OPminus3=0;OPplus3=0

for (m=0; m<5; m=m+1)
	for(t=0;t<tmax;t=t+1)
		for(i=0;i<N;i=i+1)
		phicurrent[i]=phi[i][t]-2*pi*m*(i)/N
		sinsum=sinsum+sin(phicurrent[i])
		cossum=cossum+cos(phicurrent[i])
		sinsum2=sinsum2+sin(2*phicurrent[i])
		cossum2=cossum2+cos(2*phicurrent[i])
		phicurrent[i]=phi[i][t]-(2*pi*m-pi)*i/N
		sinsum3=sinsum3+sin(2*phicurrent[i])
		cossum3=cossum3+cos(2*phicurrent[i])
		endfor
	OPminus=sqrt(sinsum^2+cossum^2)/N
	OPminus2=sqrt(abs(sqrt(sinsum2^2+cossum2^2)/N - OPminus)^2)
	OPminus3=sqrt(abs(sqrt(sinsum3^2+cossum3^2)/N - OPminus)^2)
	sinsum=0;cossum=0
	sinsum2=0;cossum2=0
	sinsum3=0;cossum3=0
		for(i=0;i<N;i=i+1)
		phicurrent[i]=phi[i][t]+2*pi*m*(i)/N
		sinsum=sinsum+sin(phicurrent[i])
		cossum=cossum+cos(phicurrent[i])
		sinsum2=sinsum2+sin(2*phicurrent[i])
		cossum2=cossum2+cos(2*phicurrent[i])
		phicurrent[i]=phi[i][t]+(2*pi*m-pi)*i/N
		sinsum3=sinsum3+sin(2*phicurrent[i])
		cossum3=cossum3+cos(2*phicurrent[i])
		endfor
	OPplus=sqrt(sinsum^2+cossum^2)/N
	OPplus2=sqrt(abs(sqrt(sinsum2^2+cossum2^2)/N - OPplus)^2)
	OPplus3=sqrt(abs(sqrt(sinsum3^2+cossum3^2)/N - OPplus)^2)
	sinsum=0;cossum=0
	sinsum2=0;cossum2=0
	sinsum3=0;cossum3=0
	OrderParameter[m][t] = max(OPplus, OPminus)
	OrderParameter2Cluster[m][t]= max(OPplus2, OPminus2)
	OrderParameter2ClusterConnected[m][t] =max(OPplus3, OPminus3)
	OPminus=0;OPplus=0
	OPminus2=0;OPplus2=0
	OPminus3=0;OPplus3=0
	endfor
endfor
end

function OPcompare(tmax)
variable tmax
variable m, currentMax
String currentString
variable OPminus, OPplus, OPminus2, OPplus2, OPminus3, OPplus3
wave  OrderParameter, OrderParameter2Cluster, OrderParameter2ClusterConnected

tmax=tmax-1
currentMax=0
for (m=0; m<5; m=m+1)

	if (currentMax<OrderParameter[m][tmax])
	currentMax=OrderParameter[m][tmax]
	currentString= "It is single cluster, mode= " + num2istr(m) +"  with Order Parameter= " + num2str(currentMax)
	endif
	
	if (currentMax<OrderParameter2Cluster[m][tmax])
	currentMax=OrderParameter2Cluster[m][tmax]
	currentString= "It is unconnected 2 cluster behavior, mode= " + num2istr(m) +"  with Order Parameter= " + num2str(currentMax)
	endif
	
	if (currentMax<OrderParameter2ClusterConnected[m][tmax] && m>0)
	currentMax=OrderParameter2ClusterConnected[m][tmax]
	currentString= "It is connected 2 cluster behavior, mode= " + num2istr(m) +"  with Order Parameter= " + num2str(currentMax)
	endif
	
endfor
Print currentString
end




