classdef Boundary < Element
    properties
        A
        B
    end
    
    methods
        function this=Boundary(A, B)
            this.A = A; this.B = B;
        end
        
        function h=plot(this) 
            h=plot([this.A(1) this.B(1)], [this.A(2) this.B(2)], 'r-'); 
        end
        
    end
    
end

