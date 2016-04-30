classdef (Abstract) Element < matlab.mixin.Heterogeneous & handle
    properties
        m=0
        V=[0 0];
    end
    
    methods (Abstract)
        plot(this)
    end
    
    methods
        function applyForce(this,F) end
        function applyForceMoment(this, R) end
        function processForces(this,dt) end
        function b=isStatic(this) 
            b=true;
        end
        function Vr=Vr(this)
            Vr=0;
        end
    end
    
end

