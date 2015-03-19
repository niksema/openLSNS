#ifndef __BIOMECH_SKETCH_H
#define __BIOMECH_SKETCH_H


//---- joint specific angles and angular velocities
theta = [x(1), x(2)-x(1)];         %[shoulder elbow]
thetaV = [x(3), x(4)-x(3)];        %[shoulder elbow]
//---- specify joint restrictions
a1min = pi/2-A1max;
a1max = pi/2-A1min;
a2min = pi-A2max;
a2max = pi-A2min;

//---- calculating total non-conservative forces
calc_ncforces();
//---- calculating total conservative forces
calc_cforces();
//---- calculating total torques in upper and lower joints;
Q1 = f1c-f1g+q1;
Q2 = f2c-f2g+q2;
//--- calculating dynamics of the system
solve_biomech( Q1, Q2 );

FB = calc_feedback( l, v, mn, f, FB );


//---- calculating total non-conservative forces
/*
%---- calculating total non-conservative forces
% 1. joints viscosity
q1v = -BJ*thetaV(1);
q2v = -BJ*thetaV(2);
% 2. joints restrictions (elastic and viscosity components)
q1e = 0;
q2e = 0;
if theta(1) > a1max && thetaV(1) > 0       % if shoulder flexion exceeds the limit
    q1e = -KJr*(theta(1)-a1max)-BJr*thetaV(1);
elseif theta(1) < a1min && thetaV(1) < 0   % if shoulder extension exceeds the limit
    q1e = -KJr*(theta(1)-a1min)-BJr*thetaV(1);
end
if theta(2) > a2max && thetaV(2) > 0        % if elbow flexion exceeds the limit
    q2e = -KJr*(theta(2)-a2max)-BJr*thetaV(2);
elseif theta(2) < a2min && thetaV(2) < 0    % if elbow flexion exceeds the limit
    q2e = -KJr*(theta(2)-a2min)-BJr*thetaV(2);
end
%3. muscles forces
% calculate motoneuronal output
X = mapFB2X( FB, X );
X = mapDr2X( DR, X );
X = mapY2X( Y, X );
Y = sp_cord( Y, B, W, X );
mn = mapY2MN( Y );
*/

//---- calculating muscle force
void calc_mforce( void )
{
    [l, v, h_up, h_dw] = geometry_block( theta, thetaV ); %get muscle geomerty
    f = muscles( mn, l, v );
    double q1m = sum(f.*h_up);
    double q2m = sum(f.*h_dw);
}

void calc_ncforces( void )
{
    //--- joints viscosity
    double q1v = -BJ*thetaV(1);
    double q2v = -BJ*thetaV(2);
    //--- joints restrictions (elastic and viscosity components)
    double q1e = 0;
    double q2e = 0;
    if theta(1) > a1max && thetaV(1) > 0       % if shoulder flexion exceeds the limit
        q1e = -KJr*(theta(1)-a1max)-BJr*thetaV(1);
    elseif theta(1) < a1min && thetaV(1) < 0   % if shoulder extension exceeds the limit
        q1e = -KJr*(theta(1)-a1min)-BJr*thetaV(1);
    end
    if theta(2) > a2max && thetaV(2) > 0        % if elbow flexion exceeds the limit
        q2e = -KJr*(theta(2)-a2max)-BJr*thetaV(2);
    elseif theta(2) < a2min && thetaV(2) < 0    % if elbow flexion exceeds the limit
        q2e = -KJr*(theta(2)-a2min)-BJr*thetaV(2);
    end
    //---- calculating muscle force
    calc_mforce();
    //---- calculating total non-conservative forces
    double q1 = q1v-q2v+q1e-q2e+q1m-q2m; %upper joint
    double q2 = q2v+q2e+q2m;             %lower joint
}

//---- calculating total conservative forces
void calc_cforces( void )
{
    double f1c = -M2*L1*D2*sin(x(1)-x(2))*x(4)^2;
    double f2c = M2*L1*D2*sin(x(1)-x(2))*x(3)^2;

    double f1g = M1*g*D1g*sin(x(1))+M2*g*L1g*sin(x(1));
    double f2g = M2*g*D2g*sin(x(2));
}

////// solve the system of dynamic equations for biomechanics
//---- input variables:
// q1, q2 - total torques in upper and lower joints, correspondingly
// I1,I2,M1,M2,L1,D2 parameters of mechanical system
// x(1), x(2) - angles between vertical line and upper/lower segments
// x(3), x(4) - derivative (angular velocity) of x(1), x(2)
//---- output variables:
// x(1), x(2) - angles between vertical line and upper/lower segments
// x(3), x(4) - derivative (angular velocity) of x(1), x(2)
void solve_biomech( double q1, double q2 )
{
    double a1 = I1+M2*L1*L1; // calculate in advance
    double a2 = I2+M2*D2*D2; // calculate in advance
    double b = M2*L1*D2*cos(x(1)-x(2));
    double D = a1*a2-b*b;

    xdot(1) = x(3); % angle in upper joint
    xdot(2) = x(4); % angle in lower joint
    xdot(3) = (a2*q1-b*q2)/D; % angular velocity in upper joint
    xdot(4) = (a1*q2-b*q1)/D; % angular velocity in lower joint
}

#endif // __BIOMECH_SKETCH_H
