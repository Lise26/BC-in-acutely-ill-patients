% network class

% Deals with networks, their construction and parameters extraction

classdef network
    properties
      Nodes {mustBeNumeric}
      Edges
      Graph
    end
    
    methods
        % Constructor
        function obj = network(matrix)
            obj.Nodes = size(matrix,1);
            obj.Edges = matrix+matrix';
            for i = 1:size(matrix,1)
                obj.Edges(i,i) = nan;
            end
        end
    end
end
