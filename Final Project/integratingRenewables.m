%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variable Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numOriginalNodes = 1700;
numAdditionalRenewableNodes = 100;     %this will be the number of new renewable energy sources we are adding to existing grid 
finalNumNodes = numOriginalNodes+numAdditionalRenewableNodes;   %total nodes when network is fully created
numEdges = 1800;
numExtraEdges = 100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End Variable Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start non-renewable energy network generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We Assume every node in the network is connected to the network in some
%way, and none of the nodes are connected to themselves (no self loops)
%The Degrees will otherwise be randomly distributed yielding a degree
%distribution that is approximately uniform 

connectionProbability = 1/numOriginalNodes;
flagForConnectionMade = 0;

A = zeros(numOriginalNodes);

%disp(A)

myVector = 1:numOriginalNodes;
%disp(myVector)

perm = myVector(randperm(length(myVector)));
%disp(perm)

for i = 1:numOriginalNodes
    while flagForConnectionMade == 0
        perm = myVector(randperm(length(myVector)));
        for j = 1:numOriginalNodes
            randomIndex = perm(j);
            if i ~= randomIndex          %doing this to avoid self loops dont want to connect to itself
                myRandNum = rand;
                if myRandNum < connectionProbability
                    A(i,randomIndex) = 1;
                    flagForConnectionMade =1;
                    break                        %break me out of inner for loop
                end                
            end 
        end  
    end
    flagForConnectionMade = 0;      %resetting flag  
end

%disp(A)
%AT THIS POINT EVERY NODE IN THE GRAPH HAS BEEN CONNECTED TO 1 OTHER RANDOM
%NODE IN THE GRAPH
extraNodeCounter = 0;
while extraNodeCounter<=numExtraEdges              %while the number of nodes is less than 100 (to be consistent with 1800 total connections)
    indices = datasample(myVector,2);               %take a uniform random sample from all of the nodes
    if indices(1) ~= indices(2)                     %if it is not about to connect to itself
        if A(indices(1),indices(2)) ~= 1            %if there is not already a connection
            A(indices(1),indices(2)) = 1;           %create a new connection between the nodes
            extraNodeCounter = extraNodeCounter+1;          
        end  
    end 
end
%disp(A)
%AT THIS POINT THE CURRENT GRID FOR CHINA HAS BEEN CREATED AS UNDIRECTED SO
%I MUST MAKE IT SYMMETRICAL AS POWER GRIDS ARE
for i = 1:numOriginalNodes
    for j = 1:numOriginalNodes
        if A(i,j) ==1
            A(j,i) =1;
        end
    end
end
%NOW THE AN UNDIRECTED GRID OF CHINA HAS BEEN CREATED
%disp(A)
firstSum = sum(A);
avgDegree = sum(firstSum)/numOriginalNodes;                 %this is <k>
%/numNodes;
disp(sum(firstSum))
disp(avgDegree)

%disp(mySparseA)

ChineseGraph = graph(A);
plot(ChineseGraph)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nonrenewable portion of network successfuly generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start adding renewables (create different strategies for connections here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end adding renewables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start finding various properties of specific network (eigenvector
%centrality etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end finding various properties of specific network (eigenvector
%centrality etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start calculating stability using our stability metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We define the stability metric of our networks as being the time averaged
%difference between the normal frequency of 60 hertz and the average
%frequency of all oscillators at that time of the network 


%WHEN A RENEWABLE ENERGY SOURCE GOES OUT ITS FREQUENCY DROPS FROM 60 to 0
%IMMEDIATELY BECAUSE THERE IS NO INERTIA

%WHEN A NON RENEWABLE ENERGY SOURCE GOES OUT ITS FREQUENCY INCREASES THE
%LONGER IT IS DISCONNECTED (MUCH SMALLER PERTURBATION) WE WILL ASSUME THAT
%IT WILL BE FIXED ON AVERAGE WHEN ITS FREQUENCY HAS INCREASED FROM 60 TO
%60.5 HERTZ

probabilityRenewableOn = 0.99;
K = 3;      %coupling constant (this will only come into play when the nodes are connected on the grid)
dampingConstant = 0.35;
numTimeSteps = 3000;           %this is how long the simulations will run
changeTime = 0.01;
stabilityTimeVector = zeros(1, numTimeSteps);      %this will hold the stability metric of the network at each timestep (averaged at the end) 
singleOscillatorFrequency = zeros(1,numTimeSteps);
avgFrequencyFrom0 = zeros(1,numTimeSteps); 
coupledComponent1Timestep = zeros(1,numTimeSteps);
avgTheta = zeros(1,numTimeSteps);
stableWn = zeros(finalNumNodes,1);     %this is what I subtract from the resulting frequencies at each timestep when finding second norm

thetaVector = zeros(finalNumNodes,1);
dthetaVector = zeros(finalNumNodes,1);

WnVector = zeros(finalNumNodes,1);              %this is wN WHICH IS WHAT WE ARE CONCERNED ABOUT
dWnVector = zeros(finalNumNodes,1);

distFromTargetWn = zeros(finalNumNodes, 1);       %this will be used to keep track of how far each generator is from target frequency

powerGeneratedVector = randn(finalNumNodes,1);      %this is the power that each oscillator/generator consumes or produces
powerGeneratedVector = powerGeneratedVector-mean(powerGeneratedVector);      %this enforces that power in=power out on grid

dampingVector = dampingConstant*ones(finalNumNodes,1);

numTrials = 1;



for trial = 1:numTrials
    if trial ==1            %here I connect the new renewables in the naive way with only 1 connection
        A(end+1:finalNumNodes,end+1:finalNumNodes) = 0;
     
        myFinalVector = 1:finalNumNodes;
        perm = myFinalVector(randperm(length(myFinalVector)));

        for i = numOriginalNodes+1:finalNumNodes   %THIS IS WHERE WE ATTACH THE NEW NODES IN DIFFERENT WAYS FIRST HERE WE ILLUSTRATE A DEAD END
            while flagForConnectionMade == 0
                perm = myFinalVector(randperm(length(myFinalVector)));
                for j = 1:finalNumNodes
                    randomIndex = perm(j);
                    if i ~= randomIndex          %doing this to avoid self loops dont want to connect to itself
                        myRandNum = rand;
                        if myRandNum < connectionProbability
                            A(i,randomIndex) = 1;
                            flagForConnectionMade =1;
                            break                        %break me out of inner for loop
                        end                
                    end 
                end  
            end
            flagForConnectionMade = 0;      %resetting flag
        end
        %MUST MAKE SYMMETRICAL AGAIN
        for i = 1:finalNumNodes
            for j = 1:finalNumNodes
                if A(i,j) ==1
                    A(j,i)=1;
                end
            end
        end        
        mySparseA = sparse(A);              %Now A will have been created with the additional renewables
    end
    
    if trial == 2
        disp("This will be running the experiment with no renewables added")
        mySparseA = sparse(A);        
    end
    
    %must setup initial conditions of network (everything stable and good
    %initially at t=0
    
    %WnVector = stableWn;                   %PUT THIS ONE BACK JACKSON
    WnVector = 0.01*randn(finalNumNodes,1);
    WnVector(1)=50;
    omega = 0.01*randn(finalNumNodes,1);
    
    
    
    for i = 1:finalNumNodes          %here I setup the initial theta values of the generators (assuming uniform dist)
        initialAngle = rand()*2*pi;
        thetaVector(i) = initialAngle;       
    end
    
    
    for t = 1:numTimeSteps
        %Must first specify what A will be for this time iteration based on
        %which renewables are on
        currA = mySparseA;          %mySparseA will stay constant and this is because generators get fixed at each timestep
        
        for i = numOriginalNodes:finalNumNodes
            myRandNum = rand;
            if abs(WnVector(i))<1          %if the frequency of the renewable generator has recovered to grid frequency
                if myRandNum>probabilityRenewableOn  %if this renewable generator is within reasonable frequency of grid and suddenly goes off
                    currA(i,:) = 0;
                    currA(:,i) = 0;
                    WnVector(i) = -60;
                end
                
            end
  
        end
        
        
        
        %%%%%%%%%%%%%%%%%%
        dthetaVector = WnVector;            %change in theta = Wn
        dWnVector = powerGeneratedVector-dampingVector.*WnVector + imag(exp(-1i*thetaVector).*(currA*exp(1i*thetaVector)));      %change in frequency
        
        %omegap = P - alpha.*omega + imag(exp(-1i*theta).*(a*exp(1i*theta)));
        
        %ABOVE I NEED TO CALCULATE dWnVector for this timestep
                      
        singleOscillatorFrequency(t) = WnVector(1550);
        
        avgFrequencyFrom0(t) = mean(WnVector.*powerGeneratedVector);
        
        
        WnVector = WnVector+changeTime*dWnVector;
        thetaVector = thetaVector+changeTime*dthetaVector;
        %disp("wankers")    
    end
    
end

plot(avgFrequencyFrom0)
%plot(coupledComponent1Timestep)
%plot(stabilityTimeVector)
%plot(singleOscillatorFrequency)
title('Frequency of a Single Nonrenewable Generator As Time Progresses')
ylabel('Frequency in Hertz')
xlabel('Time')
%plot(avgFrequencyFrom0)

%scatter(WnVector, powerGeneratedVector)

%plot stability against time
finalStabilityValue = mean(avgFrequencyFrom0);
disp("below is the final stability value for this network:")
disp(finalStabilityValue)




%STEPS TO DO:
%increase size of A to simulate adding in renewable energy sources

%run dynamics on the new matrix A where each row is an oscillator (kuramoto
%model)
%repeatedly alter A (aka cutting a connection in the grid and then test the
%stability)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





