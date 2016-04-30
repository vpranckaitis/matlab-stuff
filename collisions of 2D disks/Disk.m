classdef Disk < Element
    
    properties
        g=9.81;
        r
        P
        phi=0;
        omega=0;
        F
        M
    end
    
    methods
        function this = Disk(r, m, P, V)
            this.r=r;
            this.m=m;
            this.P=P;
            this.V=V;
            this.F=[0 -m]*this.g;
            this.M=0;
        end
        
        function h=plot(this) 
            T = [0:(pi/64):(2*pi)] + this.phi;
            X = [this.r*cos(T) 0] + this.P(1);
            Y = [this.r*sin(T) 0] + this.P(2);
            h=plot(X, Y, 'b-'); 
        end
        
        function applyForce(this, F)
            this.F = this.F + F;
        end
        
        function applyForceMoment(this, M)
            this.M = this.M + M;
        end
        
        function processForces(this, dt)
            A=this.F/this.m;
            this.V=this.V + A*dt;
            this.P=this.P + this.V * dt;
            this.F=[0 -this.m]*this.g;
            
            I=(this.m*this.r^2)/2;
            alpha=this.M/I;
            this.omega=this.omega + alpha*dt;
            this.phi=this.phi + this.omega*dt;
            this.M = 0;
        end
        
        function b=isStatic(this)
            b=false;
        end
        
        function Vr=Vr(this)
            Vr=this.r*this.omega;
        end
    end
end

