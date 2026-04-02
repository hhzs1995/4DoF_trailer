%% Nonlinear governing equations for two-wheeled trailers
% 4 DoF model, q = [psi, theta, fi, u]

% based on: H. Z. Horvath, D. Takacs, M. Antali: 
% Effects of differential braking and tire contact modeling 
% on the stability and nonlinear dynamics of towed two-wheeled trailers 
% (2026). Manuscript submitted for publication.

function dydt = trailer_4DoF_goveq(t,y,Z)

    % parameters  
    v = 60; % towing speed, m/s
    a = 0.24; % horizontal position of CoG, m
    h = 0.27; % vertical position of CoG, m
    l = 3.77; % caster length, m
    h0 = 0.35; % height of the kingpin, m
    w = 0.95; % half of the track width, m
    m = 879; % mass of the trailer, kg
    k = 30000; % stiffness of the suspension, N/m
    c = 950; % damping of the suspension, Ns/m
    klat = 13000; % lateral stiffness, N/m
    clat = 750; % lateral damping, Ns/m
    eps = 1e-3; % smoothing parameter
    Lmin = 0.28; % smallest length of the spring, m
    p = 0.065; % prestress of the spring, m (relationship: Lmax = L0-p)
	Kdpsi = 10000; % derivative gain for the braking forces, Ns
    psid0 = 0.1; % dead zone for the brake forces, rad/s
    eps_braking = 1e-3; % smoothing parameter for braking forces
    g = 9.81; % gravitational accerelation, m/s^2
    B = 32.45; C = 1.25; D = 0.65; E = -1.24; % Magic Formula parameters (constants)

    % state variables
    psi = y(1);
    theta = y(2);
    fi = y(3);
    u = y(4);
    psid = y(5);
    thetad = y(6);
    fid = y(7);
    ud = y(8);
    ylag1 = Z(:,1);
    psid_delayed = ylag1(5);

    % free length of the springs, deformations
    LB0 = (-(a*g*m)+g*l*m+2*h0*l*k)/(2*l*k);
    LJ0 = LB0;
    lambdaTB = 1/cos(fi)*(h0/cos(theta)+w*sin(fi)+l*tan(theta));
    lambdaTJ = 1/cos(fi)*(h0/cos(theta)-w*sin(fi)+l*tan(theta));

    % unit vectors (lateral, longitudinal)
    j3R0XY = [-(cos(fi)*sin(psi))+cos(psi)*sin(fi)*sin(theta);cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta);0];
    ylat = j3R0XY/norm(j3R0XY);
    xlong = [ylat(2); -ylat(1); ylat(3)];

    % kinematics
    vTBR0 = [v-w*psid*cos(fi)*cos(psi)+l*psid*cos(theta)*sin(psi)-h0*psid*cos(psi)*1/cos(theta)*tan(fi)-w*psid*cos(psi)*sin(fi)*tan(fi)-thetad*1/cos(theta)^2*(cos(psi)*(h0+l*sin(theta))+sin(psi)*(l+h0*sin(theta))*tan(fi))+h0*psid*sin(psi)*tan(theta)+l*psid*sin(psi)*sin(theta)*tan(theta)-l*psid*cos(psi)*tan(fi)*tan(theta)-fid*1/cos(fi)^2*sin(psi)*(h0*1/cos(theta)+w*sin(fi)+l*tan(theta));
        ud-l*psid*cos(psi)*cos(theta)+fid*h0*cos(psi)*1/cos(theta)+w*fid*cos(psi)*sin(fi)-w*psid*cos(fi)*sin(psi)-h0*psid*1/cos(theta)*sin(psi)*tan(fi)-w*psid*sin(fi)*sin(psi)*tan(fi)+fid*h0*cos(psi)*1/cos(theta)*tan(fi)^2+w*fid*cos(psi)*sin(fi)*tan(fi)^2-thetad*1/cos(theta)^2*(sin(psi)*(h0+l*sin(theta))-cos(psi)*(l+h0*sin(theta))*tan(fi))+fid*l*cos(psi)*tan(theta)-h0*psid*cos(psi)*tan(theta)-l*psid*cos(psi)*sin(theta)*tan(theta)-l*psid*sin(psi)*tan(fi)*tan(theta)+fid*l*cos(psi)*tan(fi)^2*tan(theta);
        0];
    vTJR0 = [v+w*psid*cos(fi)*cos(psi)+l*psid*cos(theta)*sin(psi)-h0*psid*cos(psi)*1/cos(theta)*tan(fi)+w*psid*cos(psi)*sin(fi)*tan(fi)-thetad*1/cos(theta)^2*(cos(psi)*(h0+l*sin(theta))+sin(psi)*(l+h0*sin(theta))*tan(fi))+h0*psid*sin(psi)*tan(theta)+l*psid*sin(psi)*sin(theta)*tan(theta)-l*psid*cos(psi)*tan(fi)*tan(theta)-fid*1/cos(fi)^2*sin(psi)*(h0*1/cos(theta)-w*sin(fi)+l*tan(theta));
        ud-l*psid*cos(psi)*cos(theta)+fid*h0*cos(psi)*1/cos(theta)-w*fid*cos(psi)*sin(fi)+w*psid*cos(fi)*sin(psi)-h0*psid*1/cos(theta)*sin(psi)*tan(fi)+w*psid*sin(fi)*sin(psi)*tan(fi)+fid*h0*cos(psi)*1/cos(theta)*tan(fi)^2-w*fid*cos(psi)*sin(fi)*tan(fi)^2-thetad*1/cos(theta)^2*(sin(psi)*(h0+l*sin(theta))-cos(psi)*(l+h0*sin(theta))*tan(fi))+fid*l*cos(psi)*tan(theta)-h0*psid*cos(psi)*tan(theta)-l*psid*cos(psi)*sin(theta)*tan(theta)-l*psid*sin(psi)*tan(fi)*tan(theta)+fid*l*cos(psi)*tan(fi)^2*tan(theta);
        0];

    vTBR3 = [(v*cos(psi)*cos(theta)-w*psid*cos(theta)*1/cos(fi)+v*cos(psi)*1/cos(theta)-w*psid*1/cos(fi)*1/cos(theta)+2*ud*cos(theta)*sin(psi)-2*thetad*1/cos(theta)*(h0+l*sin(theta))-2*h0*psid*tan(fi)-2*l*psid*sin(theta)*tan(fi)-v*cos(psi)*sin(theta)*tan(theta)+w*psid*1/cos(fi)*sin(theta)*tan(theta))/2;
        ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+(-(l*psid*cos(fi)*cos(theta))+l*psid*cos(theta)*1/cos(fi)-l*psid*cos(fi)*1/cos(theta)-3*l*psid*1/cos(fi)*1/cos(theta)+4*l*thetad*sin(fi)-2*v*cos(fi)*sin(psi)-2*v*1/cos(fi)*sin(psi)+4*v*cos(psi)*sin(fi)*sin(theta)+l*psid*cos(theta)*sin(fi)*tan(fi)+l*psid*1/cos(theta)*sin(fi)*tan(fi)+2*v*sin(fi)*sin(psi)*tan(fi)-4*w*psid*sin(theta)*tan(fi)-4*h0*psid*1/cos(fi)*tan(theta)+l*psid*cos(fi)*sin(theta)*tan(theta)-l*psid*1/cos(fi)*sin(theta)*tan(theta)-l*psid*sin(fi)*sin(theta)*tan(fi)*tan(theta)+4*fid*(h0*1/cos(fi)*1/cos(theta)+w*tan(fi)+l*1/cos(fi)*tan(theta)))/4;
        ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+(thetad*1/cos(fi)*1/cos(theta)^2*(-6*l+2*l*cos(2*fi)+2*l*cos(2*theta)+l*cos(2*(-fi+theta))+l*cos(2*(fi+theta))-8*h0*sin(theta))+8*(l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)-w*psid*sin(theta)+v*cos(fi)*cos(psi)*sin(theta)-fid*1/cos(fi)*tan(fi)*(h0*1/cos(theta)+w*sin(fi)+l*tan(theta))))/8];
    vTJR3 = [(v*cos(psi)*cos(theta)+w*psid*cos(theta)*1/cos(fi)+v*cos(psi)*1/cos(theta)+w*psid*1/cos(fi)*1/cos(theta)+2*ud*cos(theta)*sin(psi)-2*thetad*1/cos(theta)*(h0+l*sin(theta))-2*h0*psid*tan(fi)-2*l*psid*sin(theta)*tan(fi)-v*cos(psi)*sin(theta)*tan(theta)-w*psid*1/cos(fi)*sin(theta)*tan(theta))/2;
        ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+(-(l*psid*cos(fi)*cos(theta))+l*psid*cos(theta)*1/cos(fi)-l*psid*cos(fi)*1/cos(theta)-3*l*psid*1/cos(fi)*1/cos(theta)+4*l*thetad*sin(fi)-2*v*cos(fi)*sin(psi)-2*v*1/cos(fi)*sin(psi)+4*v*cos(psi)*sin(fi)*sin(theta)+l*psid*cos(theta)*sin(fi)*tan(fi)+l*psid*1/cos(theta)*sin(fi)*tan(fi)+2*v*sin(fi)*sin(psi)*tan(fi)+4*w*psid*sin(theta)*tan(fi)-4*h0*psid*1/cos(fi)*tan(theta)+l*psid*cos(fi)*sin(theta)*tan(theta)-l*psid*1/cos(fi)*sin(theta)*tan(theta)-l*psid*sin(fi)*sin(theta)*tan(fi)*tan(theta)+4*fid*(h0*1/cos(fi)*1/cos(theta)-w*tan(fi)+l*1/cos(fi)*tan(theta)))/4;
        ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+(thetad*1/cos(fi)*1/cos(theta)^2*(-6*l+2*l*cos(2*fi)+2*l*cos(2*theta)+l*cos(2*(-fi+theta))+l*cos(2*(fi+theta))-8*h0*sin(theta))+8*(l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+w*psid*sin(theta)+v*cos(fi)*cos(psi)*sin(theta)-fid*1/cos(fi)*tan(fi)*(h0*1/cos(theta)-w*sin(fi)+l*tan(theta))))/8];

    vBR3 = [(-(w*psid*cos(fi))+v*cos(psi))*cos(theta)+w*thetad*sin(fi)+ud*cos(theta)*sin(psi);
        -(l*psid*cos(fi)*cos(theta))+l*thetad*sin(fi)-v*cos(fi)*sin(psi)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta));
        w*fid+l*thetad*cos(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)-w*psid*sin(theta)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))];
    vJR3 = [(w*psid*cos(fi)+v*cos(psi))*cos(theta)-w*thetad*sin(fi)+ud*cos(theta)*sin(psi);
        -(l*psid*cos(fi)*cos(theta))+l*thetad*sin(fi)-v*cos(fi)*sin(psi)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta));
        -(w*fid)+l*thetad*cos(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+w*psid*sin(theta)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))];

    % side slip angles
    alfaB = -atan(vTBR0'*ylat/sqrt((vTBR0-(vTBR0'*ylat)*ylat)'*(vTBR0-(vTBR0'*ylat)*ylat)));
    alfaJ = -atan(vTJR0'*ylat/sqrt((vTJR0-(vTJR0'*ylat)*ylat)'*(vTJR0-(vTJR0'*ylat)*ylat)));

    % forces in the suspension system
    kInf = 2*k;
    cInf= 2*c;
    FBrugo10 = (LB0-lambdaTB)*k+(vTBR3(3)-vBR3(3))*c;
    FBrugo20 = (Lmin-lambdaTB)*kInf+(vTBR3(3)-vBR3(3))*cInf;	
	FBrugo = (FBrugo10*(1+tanh((LB0-p-lambdaTB)/eps))/2+FBrugo20*(1+tanh((Lmin-lambdaTB)/eps))/2)*(1+tanh((FBrugo10*(1+tanh((LB0-p-lambdaTB)/eps))/2+FBrugo20*(1+tanh((Lmin-lambdaTB)/eps))/2)/eps))/2;
    FJrugo10 = (LJ0-lambdaTJ)*k+(vTJR3(3)-vJR3(3))*c;
    FJrugo20 = (Lmin-lambdaTJ)*kInf+(vTJR3(3)-vJR3(3))*cInf;
	FJrugo = (FJrugo10*(1+tanh((LJ0-p-lambdaTJ)/eps))/2+FJrugo20*(1+tanh((Lmin-lambdaTJ)/eps))/2)*(1+tanh((FJrugo10*(1+tanh((LJ0-p-lambdaTJ)/eps))/2+FJrugo20*(1+tanh((Lmin-lambdaTJ)/eps))/2)/eps))/2;

    % Magic Formula
    MFJ = D*sin(C*atan(B*alfaJ-E*(B*alfaJ-atan(B*alfaJ))));
    MFB = D*sin(C*atan(B*alfaB-E*(B*alfaB-atan(B*alfaB))));

    % brake forces with Heaviside-functions, considering the dead zone of the controller
    FBbrake = -Kdpsi*(psid_delayed+psid0)*(1-tanh((psid_delayed+psid0)/eps_braking))/2; 
    FJbrake = Kdpsi*(psid_delayed-psid0)*(1+tanh((psid_delayed-psid0)/eps_braking))/2;
    FJbrake0 = -FJbrake*xlong;
    FBbrake0 = -FBbrake*xlong;
    FJbrake0x = FJbrake0(1);
    FJbrake0y = FJbrake0(2);
    FJbrake0z = FJbrake0(3);
    FBbrake0x = FBbrake0(1);
    FBbrake0y = FBbrake0(2);
    FBbrake0z = FBbrake0(3);
    
    % normal forces
    FTJZ0 = (1/cos(theta)*sqrt(cos(fi)^2+sin(fi)^2*sin(theta)^2)*(-FJbrake0z+1/cos(theta)*(FJrugo*1/cos(fi)+(FJbrake0y*cos(psi)-FJbrake0x*sin(psi))*tan(fi))-(FJbrake0x*cos(psi)+FJbrake0y*sin(psi))*tan(theta)))/(-(MFJ*sin(fi))+1/cos(theta)*sqrt(cos(fi)^2+sin(fi)^2*sin(theta)^2));
    FTBZ0 = (1/cos(theta)*sqrt(cos(fi)^2+sin(fi)^2*sin(theta)^2)*(-FBbrake0z+1/cos(theta)*(FBrugo*1/cos(fi)+(FBbrake0y*cos(psi)-FBbrake0x*sin(psi))*tan(fi))-(FBbrake0x*cos(psi)+FBbrake0y*sin(psi))*tan(theta)))/(-(MFB*sin(fi))+1/cos(theta)*sqrt(cos(fi)^2+sin(fi)^2*sin(theta)^2));

    % lateral tyre forces
    FJtyre = FTJZ0*MFJ;
    FBtyre = FTBZ0*MFB;
    FJtyre0 = FJtyre*ylat;
    FBtyre0 = FBtyre*ylat;
    FJtyre0x = FJtyre0(1);
    FJtyre0y = FJtyre0(2);
    FJtyre0z = FJtyre0(3);
    FBtyre0x = FBtyre0(1);
    FBtyre0y = FBtyre0(2);
    FBtyre0z = FBtyre0(3);

    % force at the kingpin
    FArugo = -klat*u-clat*ud;

    % mass moment of intertia
    Thetax = ((4*w^2+4*h^2)*m)/6; 
    Thetay = ((4*h^2+l^2)*m)/6; 
    Thetaz = ((4*w^2+l^2)*m)/6;

    % mass matrix
    M11 = Thetaz*cos(fi)^2*cos(theta)^2+h^2*m*cos(theta)^2*sin(fi)^2+Thetay*cos(theta)^2*sin(fi)^2-a*m*cos(theta)*sin(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))+l*m*cos(theta)*sin(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))+Thetax*sin(theta)^2+a*m*cos(fi)*cos(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))-l*m*cos(fi)*cos(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))+h*m*sin(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta));
    M12 = h^2*m*cos(fi)*cos(theta)*sin(fi)+Thetay*cos(fi)*cos(theta)*sin(fi)-Thetaz*cos(fi)*cos(theta)*sin(fi)-(a-l)*m*cos(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))-(a-l)*m*sin(fi)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta));
    M13 = -(Thetax*sin(theta))-h*m*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta));
    M14 = h*m*cos(theta)^2*sin(fi)*sin(psi)+m*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+m*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta));
    M21 = h^2*m*cos(fi)*cos(theta)*sin(fi)+Thetay*cos(fi)*cos(theta)*sin(fi)-Thetaz*cos(fi)*cos(theta)*sin(fi)-h*(a-l)*m*sin(fi)*sin(theta);
    M22 = h^2*m*cos(fi)^2+(a-l)^2*m*cos(fi)^2+Thetay*cos(fi)^2+(a-l)^2*m*sin(fi)^2+Thetaz*sin(fi)^2;
    M23 = h*(a-l)*m*sin(fi);
    M24 = h*m*cos(fi)*cos(theta)*sin(psi)-(a-l)*m*cos(fi)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))-(a-l)*m*sin(fi)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta));
    M31 = -(a*h*m*cos(fi)*cos(theta))+h*l*m*cos(fi)*cos(theta)-h^2*m*sin(theta)-Thetax*sin(theta);
    M32 = M23;
    M33 = h^2*m+Thetax;
    M34 = -(h*m*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta)));
    M41 = h*m*cos(theta)^2*sin(fi)*sin(psi)-a*m*cos(theta)*sin(fi)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+l*m*cos(theta)*sin(fi)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+a*m*cos(fi)*cos(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))-l*m*cos(fi)*cos(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+h*m*sin(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta));
    M42 = M24;
    M43 = M34;
    M44 = m*cos(theta)^2*sin(psi)^2+m*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))^2+m*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))^2;

    MM = [M11,M12,M13,M14; 
          M21,M22,M23,M24;
          M31,M32,M33,M34;
          M41,M42,M43,M44];

    % vector containing nonlinearities
    Cv1 = -(fid*thetad*Thetaz*cos(fi)^2*cos(theta))+fid*h^2*m*psid*cos(fi)*cos(theta)^2*sin(fi)+fid*psid*Thetay*cos(fi)*cos(theta)^2*sin(fi)-fid*psid*Thetaz*cos(fi)*cos(theta)^2*sin(fi)+h*m*psid*ud*cos(psi)*cos(theta)^2*sin(fi)-fid*h^2*m*thetad*cos(theta)*sin(fi)^2-fid*thetad*Thetay*cos(theta)*sin(fi)^2+fid*Thetaz*cos(theta)*sin(fi)*(-(psid*cos(fi)*cos(theta))+thetad*sin(fi))-a*fid*m*psid*cos(fi)*cos(theta)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))+fid*l*m*psid*cos(fi)*cos(theta)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))+m*thetad*v*cos(fi)*cos(psi)*cos(theta)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))+fid*(a-l)*m*thetad*sin(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))+m*psid*v*cos(psi)*sin(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))+fid*Thetay*cos(fi)*cos(theta)*(thetad*cos(fi)+psid*cos(theta)*sin(fi))-h*m*psid*v*cos(theta)^2*sin(fi)*sin(psi)+fid*m*v*cos(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))*sin(psi)+fid*h*m*cos(fi)*cos(theta)*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))+psid*thetad*Thetax*cos(theta)*sin(theta)-psid*thetad*Thetaz*cos(fi)^2*cos(theta)*sin(theta)-psid*thetad*Thetay*cos(theta)*sin(fi)^2*sin(theta)-h*m*thetad*cos(theta)*sin(fi)*(v*cos(psi)+h*psid*sin(fi))*sin(theta)+thetad*Thetaz*cos(fi)*(-(psid*cos(fi)*cos(theta))+thetad*sin(fi))*sin(theta)+a*m*psid*thetad*sin(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))*sin(theta)-l*m*psid*thetad*sin(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))*sin(theta)-fid*m*v*cos(psi)*sin(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))*sin(theta)-thetad*Thetay*sin(fi)*(thetad*cos(fi)+psid*cos(theta)*sin(fi))*sin(theta)-h*m*thetad*ud*cos(theta)*sin(fi)*sin(psi)*sin(theta)-m*psid*v*cos(fi)*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))*sin(psi)*sin(theta)-h*m*thetad*sin(fi)*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))*sin(theta)-fid*(a-l)*m*thetad*cos(fi)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))-m*psid*v*cos(fi)*cos(psi)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))+h*m*psid*thetad*cos(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))-a*fid*m*psid*cos(theta)*sin(fi)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))+fid*l*m*psid*cos(theta)*sin(fi)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))+m*thetad*v*cos(psi)*cos(theta)*sin(fi)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))+fid*m*v*sin(fi)*sin(psi)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))-a*m*psid*thetad*cos(fi)*sin(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))+l*m*psid*thetad*cos(fi)*sin(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))+fid*m*v*cos(fi)*cos(psi)*sin(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))-m*psid*v*sin(fi)*sin(psi)*sin(theta)*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))-thetad*Thetax*cos(theta)*(fid-psid*sin(theta))+m*ud*(a*cos(fi)*cos(theta)-l*cos(fi)*cos(theta)+h*sin(theta))*(-(fid*cos(psi)*sin(fi))-psid*cos(fi)*sin(psi)+thetad*cos(theta)*sin(fi)*sin(psi)+psid*cos(psi)*sin(fi)*sin(theta)+fid*cos(fi)*sin(psi)*sin(theta))+m*ud*(-(a*cos(theta)*sin(fi))+l*cos(theta)*sin(fi))*(-(fid*cos(fi)*cos(psi))+thetad*cos(fi)*cos(theta)*sin(psi)+psid*sin(fi)*sin(psi)+psid*cos(fi)*cos(psi)*sin(theta)-fid*sin(fi)*sin(psi)*sin(theta))+m*(-(a*fid*cos(fi)*cos(theta))+fid*l*cos(fi)*cos(theta)+a*thetad*sin(fi)*sin(theta)-l*thetad*sin(fi)*sin(theta))*(-((a-l)*thetad*cos(fi))-a*psid*cos(theta)*sin(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta)))+m*(h*thetad*cos(theta)-a*fid*cos(theta)*sin(fi)+fid*l*cos(theta)*sin(fi)-a*thetad*cos(fi)*sin(theta)+l*thetad*cos(fi)*sin(theta))*(-(fid*h)+a*psid*cos(fi)*cos(theta)-l*psid*cos(fi)*cos(theta)-(a-l)*thetad*sin(fi)-v*cos(fi)*sin(psi)+h*psid*sin(theta)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta)))-(m*(2*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))*(ud*cos(psi)*cos(theta)-v*cos(theta)*sin(psi))+2*(v*cos(psi)*sin(fi)-v*cos(fi)*sin(psi)*sin(theta)+ud*(sin(fi)*sin(psi)+cos(fi)*cos(psi)*sin(theta)))*(-((a-l)*thetad*cos(fi))-a*psid*cos(theta)*sin(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta)))+2*(-(v*cos(fi)*cos(psi))-v*sin(fi)*sin(psi)*sin(theta)+ud*(-(cos(fi)*sin(psi))+cos(psi)*sin(fi)*sin(theta)))*(-(fid*h)+a*psid*cos(fi)*cos(theta)-l*psid*cos(fi)*cos(theta)-(a-l)*thetad*sin(fi)-v*cos(fi)*sin(psi)+h*psid*sin(theta)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta)))))/2+cos(psi)*(w*(FBbrake0x+FBtyre0x-FJbrake0x-FJtyre0x)*1/cos(fi)+1/cos(theta)*((FBbrake0y+FBtyre0y+FJbrake0y+FJtyre0y)*l+(FBbrake0x+FBtyre0x+FJbrake0x+FJtyre0x)*h0*tan(fi))+((FBbrake0y+FBtyre0y+FJbrake0y+FJtyre0y)*h0+(FBbrake0x+FBtyre0x+FJbrake0x+FJtyre0x)*l*tan(fi))*tan(theta))+sin(psi)*(w*(FBbrake0y+FBtyre0y-FJbrake0y-FJtyre0y)*1/cos(fi)+1/cos(theta)*(-((FBbrake0x+FBtyre0x+FJbrake0x+FJtyre0x)*l)+(FBbrake0y+FBtyre0y+FJbrake0y+FJtyre0y)*h0*tan(fi))+(-((FBbrake0x+FBtyre0x+FJbrake0x+FJtyre0x)*h0)+(FBbrake0y+FBtyre0y+FJbrake0y+FJtyre0y)*l*tan(fi))*tan(theta));
    Cv2 = fid*h^2*m*psid*cos(fi)^2*cos(theta)+a*fid*(a-l)*m*psid*cos(fi)^2*cos(theta)-fid*(a-l)*l*m*psid*cos(fi)^2*cos(theta)+fid*psid*Thetay*cos(fi)^2*cos(theta)+h*m*psid*ud*cos(fi)*cos(psi)*cos(theta)-(a-l)*m*thetad*v*cos(fi)^2*cos(psi)*cos(theta)-fid*h^2*m*thetad*cos(fi)*sin(fi)-fid*thetad*Thetay*cos(fi)*sin(fi)+fid*thetad*Thetaz*cos(fi)*sin(fi)-h*(a-l)*m*psid*thetad*cos(theta)*sin(fi)+a*fid*(a-l)*m*psid*cos(theta)*sin(fi)^2-fid*(a-l)*l*m*psid*cos(theta)*sin(fi)^2+fid*psid*Thetaz*cos(theta)*sin(fi)^2-(a-l)*m*thetad*v*cos(psi)*cos(theta)*sin(fi)^2+fid*Thetaz*cos(fi)*(-(psid*cos(fi)*cos(theta))+thetad*sin(fi))-fid*Thetay*sin(fi)*(thetad*cos(fi)+psid*cos(theta)*sin(fi))-fid*(a-l)*m*v*cos(fi)^2*sin(psi)-h*m*psid*v*cos(fi)*cos(theta)*sin(psi)-fid*(a-l)*m*v*sin(fi)^2*sin(psi)-fid*h*m*sin(fi)*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))-psid*thetad*Thetay*cos(fi)*sin(fi)*sin(theta)+psid*thetad*Thetaz*cos(fi)*sin(fi)*sin(theta)-h*m*thetad*cos(fi)*(v*cos(psi)+h*psid*sin(fi))*sin(theta)-h*m*thetad*ud*cos(fi)*sin(psi)*sin(theta)+(a-l)*m*psid*v*cos(fi)^2*sin(psi)*sin(theta)+(a-l)*m*psid*v*sin(fi)^2*sin(psi)*sin(theta)-(a-l)*m*ud*sin(fi)*(-(fid*cos(psi)*sin(fi))-psid*cos(fi)*sin(psi)+thetad*cos(theta)*sin(fi)*sin(psi)+psid*cos(psi)*sin(fi)*sin(theta)+fid*cos(fi)*sin(psi)*sin(theta))-(a-l)*m*ud*cos(fi)*(-(fid*cos(fi)*cos(psi))+thetad*cos(fi)*cos(theta)*sin(psi)+psid*sin(fi)*sin(psi)+psid*cos(fi)*cos(psi)*sin(theta)-fid*sin(fi)*sin(psi)*sin(theta))+fid*(a-l)*m*sin(fi)*(-((a-l)*thetad*cos(fi))-a*psid*cos(theta)*sin(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta)))-fid*(a-l)*m*cos(fi)*(-(fid*h)+a*psid*cos(fi)*cos(theta)-l*psid*cos(fi)*cos(theta)-(a-l)*thetad*sin(fi)-v*cos(fi)*sin(psi)+h*psid*sin(theta)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta)))+(-4*g*(a-l)*m*cos(theta)-3*h*g*m*cos(fi)*sin(theta)-1/cos(fi)*1/cos(theta)^2*(3*(FBrugo+FJrugo)+4*(FBbrake0y+FBtyre0y+FJbrake0y+FJtyre0y)*cos(psi)*sin(fi)+(FBrugo+FJrugo)*sin(fi)^2-4*(FBbrake0x+FBtyre0x+FJbrake0x+FJtyre0x)*sin(fi)*sin(psi))*(l+h0*sin(theta))+4*1/cos(theta)^2*((FBbrake0x+FBtyre0x+FJbrake0x+FJtyre0x)*cos(psi)+(FBbrake0y+FBtyre0y+FJbrake0y+FJtyre0y)*sin(psi))*(h0+l*sin(theta))-cos(fi)*1/cos(theta)^2*((FBrugo+FJrugo)*l+(FBrugo*h0+FJrugo*h0+h*g*m)*sin(theta)-h*g*m*sin(theta)^3))/4+(-2*psid*Thetaz*cos(fi)*(-(psid*cos(fi)*cos(theta))+thetad*sin(fi))*sin(theta)+2*psid*Thetay*sin(fi)*(thetad*cos(fi)+psid*cos(theta)*sin(fi))*sin(theta)+2*psid*Thetax*cos(theta)*(fid-psid*sin(theta))-m*(2*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))*(-((v*cos(psi)+h*psid*sin(fi))*sin(theta))-ud*sin(psi)*sin(theta))+2*(v*cos(fi)*cos(psi)*cos(theta)+ud*cos(fi)*cos(theta)*sin(psi)+a*psid*sin(fi)*sin(theta)-l*psid*sin(fi)*sin(theta))*(-((a-l)*thetad*cos(fi))-a*psid*cos(theta)*sin(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta)))+2*(h*psid*cos(theta)+v*cos(psi)*cos(theta)*sin(fi)+ud*cos(theta)*sin(fi)*sin(psi)-a*psid*cos(fi)*sin(theta)+l*psid*cos(fi)*sin(theta))*(-(fid*h)+a*psid*cos(fi)*cos(theta)-l*psid*cos(fi)*cos(theta)-(a-l)*thetad*sin(fi)-v*cos(fi)*sin(psi)+h*psid*sin(theta)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta)))))/2;
    Cv3 = fid*h*(a-l)*m*thetad*cos(fi)+h*m*psid*v*cos(fi)*cos(psi)-h^2*m*psid*thetad*cos(theta)-psid*thetad*Thetax*cos(theta)+a*fid*h*m*psid*cos(theta)*sin(fi)-fid*h*l*m*psid*cos(theta)*sin(fi)-h*m*thetad*v*cos(psi)*cos(theta)*sin(fi)-fid*h*m*v*sin(fi)*sin(psi)+a*h*m*psid*thetad*cos(fi)*sin(theta)-h*l*m*psid*thetad*cos(fi)*sin(theta)-fid*h*m*v*cos(fi)*cos(psi)*sin(theta)+h*m*psid*v*sin(fi)*sin(psi)*sin(theta)-h*m*ud*(-(fid*cos(psi)*sin(fi))-psid*cos(fi)*sin(psi)+thetad*cos(theta)*sin(fi)*sin(psi)+psid*cos(psi)*sin(fi)*sin(theta)+fid*cos(fi)*sin(psi)*sin(theta))+(-2*Thetay*(psid*cos(fi)*cos(theta)-thetad*sin(fi))*(thetad*cos(fi)+psid*cos(theta)*sin(fi))-2*Thetaz*(-(psid*cos(fi)*cos(theta))+thetad*sin(fi))*(thetad*cos(fi)+psid*cos(theta)*sin(fi))-m*(2*(h*psid*cos(fi)*cos(theta)-h*thetad*sin(fi))*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))+2*(-((a-l)*thetad*cos(fi))-a*psid*cos(theta)*sin(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta)))*(-(a*psid*cos(fi)*cos(theta))+l*psid*cos(fi)*cos(theta)+(a-l)*thetad*sin(fi)+v*cos(fi)*sin(psi)-v*cos(psi)*sin(fi)*sin(theta)+ud*(-(cos(fi)*cos(psi))-sin(fi)*sin(psi)*sin(theta)))+2*(-((a-l)*thetad*cos(fi))-a*psid*cos(theta)*sin(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta)))*(-(fid*h)+a*psid*cos(fi)*cos(theta)-l*psid*cos(fi)*cos(theta)-(a-l)*thetad*sin(fi)-v*cos(fi)*sin(psi)+h*psid*sin(theta)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta)))))/2+(2*h*g*m*cos(theta)*sin(fi)*(-3+tan(fi)^2)-8*1/cos(fi)^2*(w*FBrugo-w*FJrugo-(FBbrake0x+FBtyre0x+FJbrake0x+FJtyre0x)*h0*1/cos(theta)*sin(psi)+(FBbrake0y+FBtyre0y+FJbrake0y+FJtyre0y)*cos(psi)*1/cos(theta)*(h0+l*sin(theta))-FBbrake0x*l*sin(psi)*tan(theta)-FBtyre0x*l*sin(psi)*tan(theta)-FJbrake0x*l*sin(psi)*tan(theta)-FJtyre0x*l*sin(psi)*tan(theta))-1/cos(fi)*tan(fi)*(8*w*(FBbrake0y+FBtyre0y-FJbrake0y-FJtyre0y)*cos(psi)+h*g*m*cos(theta)+8*FBrugo*h0*1/cos(theta)+8*FJrugo*h0*1/cos(theta)+h*g*m*1/cos(theta)-8*w*FBbrake0x*sin(psi)-8*w*FBtyre0x*sin(psi)+8*w*FJbrake0x*sin(psi)+8*w*FJtyre0x*sin(psi)+8*FBrugo*l*tan(theta)+8*FJrugo*l*tan(theta)-h*g*m*sin(theta)*tan(theta)))/8;
    Cv4 = -FArugo-FBbrake0y-FBtyre0y-FJbrake0y-FJtyre0y+fid*h*m*psid*cos(fi)*cos(theta)^2*sin(psi)+m*psid*ud*cos(psi)*cos(theta)^2*sin(psi)-fid*h*m*thetad*cos(theta)*sin(fi)*sin(psi)-m*psid*v*cos(theta)^2*sin(psi)^2+m*psid*cos(psi)*cos(theta)*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))-m*thetad*cos(theta)*(v*cos(psi)+h*psid*sin(fi))*sin(psi)*sin(theta)-m*thetad*ud*cos(theta)*sin(psi)^2*sin(theta)-m*thetad*sin(psi)*(h*thetad*cos(fi)+cos(theta)*(v*cos(psi)+h*psid*sin(fi))+ud*cos(theta)*sin(psi))*sin(theta)-a*fid*m*psid*cos(fi)*cos(theta)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+fid*l*m*psid*cos(fi)*cos(theta)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+m*thetad*v*cos(fi)*cos(psi)*cos(theta)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+fid*(a-l)*m*thetad*sin(fi)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+m*psid*v*cos(psi)*sin(fi)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+fid*m*v*cos(fi)*sin(psi)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))+a*m*psid*thetad*sin(fi)*sin(theta)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))-l*m*psid*thetad*sin(fi)*sin(theta)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))-fid*m*v*cos(psi)*sin(fi)*sin(theta)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))-m*psid*v*cos(fi)*sin(psi)*sin(theta)*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))-fid*(a-l)*m*thetad*cos(fi)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))-m*psid*v*cos(fi)*cos(psi)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+h*m*psid*thetad*cos(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))-a*fid*m*psid*cos(theta)*sin(fi)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+fid*l*m*psid*cos(theta)*sin(fi)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+m*thetad*v*cos(psi)*cos(theta)*sin(fi)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+fid*m*v*sin(fi)*sin(psi)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))-a*m*psid*thetad*cos(fi)*sin(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+l*m*psid*thetad*cos(fi)*sin(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+fid*m*v*cos(fi)*cos(psi)*sin(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))-m*psid*v*sin(fi)*sin(psi)*sin(theta)*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+m*ud*(-(fid*cos(psi)*sin(fi))-psid*cos(fi)*sin(psi)+thetad*cos(theta)*sin(fi)*sin(psi)+psid*cos(psi)*sin(fi)*sin(theta)+fid*cos(fi)*sin(psi)*sin(theta))*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta))+m*ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta))*(-(fid*cos(fi)*cos(psi))+thetad*cos(fi)*cos(theta)*sin(psi)+psid*sin(fi)*sin(psi)+psid*cos(fi)*cos(psi)*sin(theta)-fid*sin(fi)*sin(psi)*sin(theta))+m*(-(fid*cos(fi)*cos(psi))+thetad*cos(fi)*cos(theta)*sin(psi)+psid*sin(fi)*sin(psi)+psid*cos(fi)*cos(psi)*sin(theta)-fid*sin(fi)*sin(psi)*sin(theta))*(-((a-l)*thetad*cos(fi))-a*psid*cos(theta)*sin(fi)+l*psid*cos(theta)*sin(fi)+v*sin(fi)*sin(psi)+v*cos(fi)*cos(psi)*sin(theta)+ud*(-(cos(psi)*sin(fi))+cos(fi)*sin(psi)*sin(theta)))+m*(-(fid*cos(psi)*sin(fi))-psid*cos(fi)*sin(psi)+thetad*cos(theta)*sin(fi)*sin(psi)+psid*cos(psi)*sin(fi)*sin(theta)+fid*cos(fi)*sin(psi)*sin(theta))*(-(fid*h)+a*psid*cos(fi)*cos(theta)-l*psid*cos(fi)*cos(theta)-(a-l)*thetad*sin(fi)-v*cos(fi)*sin(psi)+h*psid*sin(theta)+v*cos(psi)*sin(fi)*sin(theta)+ud*(cos(fi)*cos(psi)+sin(fi)*sin(psi)*sin(theta)));

    Cv = [Cv1;Cv2;Cv3;Cv4];

    % nonlinear governing equations
    dydt = [psid;thetad;fid;ud;-inv(MM)*Cv];

end