numOriginalNodes = 5;
numAdditionalRenewableNodes = 2;     %this will be the number of new renewable energy sources we are adding to existing grid 
finalNumNodes = numOriginalNodes+numAdditionalRenewableNodes;   %total nodes when network is fully created
numEdges = 5;
numExtraEdges = 2;



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
for i = 1:numOriginalNodes
    for j = 1:numOriginalNodes
        if A(i,j) ==1
            A(j,i)=1;
        end
    end
end

ChineseGraph = graph(A);
plot(ChineseGraph)

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
            A(j,i)=1;
        end
    end
end

%disp(A)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A(end+1:finalNumNodes,end+1:finalNumNodes) = 0;
%disp(A)


myFinalVector = 1:finalNumNodes;
perm = myFinalVector(randperm(length(myFinalVector)));

for i = numOriginalNodes+1:finalNumNodes   %THIS IS WHERE WE ATTACH THE NEW NODES IN DIFFERENT WAYS FIRST HERE WE ILLUSTRATE A DEAD END
    connectionCount = 0;
    while connectionCount < 5
        perm = myFinalVector(randperm(length(myFinalVector)));
        for j = 1:finalNumNodes
            randomIndex = perm(j);
            if i ~= randomIndex          %doing this to avoid self loops dont want to connect to itself
                if A(i,randomIndex)~=1
                    myRandNum = rand;
                    if myRandNum < connectionProbability
                        A(i,randomIndex) = 1;
                        connectionCount =connectionCount+1;
                        break                        %break me out of inner for loop
                    end
                end                
            end 
        end  
    end
end
disp("wankers")


%disp(A)

%MUST MAKE SYMMETRICAL AGAIN
for i = 1:finalNumNodes
    for j = 1:finalNumNodes
        if A(i,j) ==1
            A(j,i)=1;
        end
    end
end
%disp(A)
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cutting a node

A(6,:) = 0;
%disp(A);

A(:,6) = 0;         %Above I cut out a node
%disp(A);


%PROCEDURE FOR DETACHING A NODE DURING A TIMESTEP--> ZERO OUT THE NODES ROW
%AND COLUMN, AND SET ITS CORRESPONDING FREQUENCY TO -60







