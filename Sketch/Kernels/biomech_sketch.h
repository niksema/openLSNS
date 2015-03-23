#ifndef __BIOMECH_SKETCH_H
#define __BIOMECH_SKETCH_H



//---- calculating total non-conservative forces
/*
    % joint specific angles and angular velocities
    theta = [x(1), x(2)-x(1)];         %[shoulder elbow]
    thetaV = [x(3), x(4)-x(3)];        %[shoulder elbow]
    % specify joint restrictions
    a1min = pi/2-A1max;
    a1max = pi/2-A1min;
    a2min = pi-A2max;
    a2max = pi-A2min;

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
    % calculate muscle force
    [l, v, h_up, h_dw] = geometry_block( theta, thetaV ); %get muscle geomerty
    f = muscles( mn, l, v );
    FB = calc_feedback( l, v, mn, f, FB );
    q1m = sum(f.*h_up);
    q2m = sum(f.*h_dw);
    %4. total non-conservative forces
    q1 = q1v-q2v+q1e-q2e+q1m-q2m; %upper joint
    q2 = q2v+q2e+q2m;             %lower joint
    %---- calculating total conservative forces
    f1c = -M2*L1*D2*sin(x(1)-x(2))*x(4)^2;
    f2c = M2*L1*D2*sin(x(1)-x(2))*x(3)^2;
    f1g = M1*g*D1g*sin(x(1))+M2*g*L1g*sin(x(1));
    f2g = M2*g*D2g*sin(x(2));
    %---- calculating total torques in upper and lower joints;
    Q1 = f1c-f1g+q1;
    Q2 = f2c-f2g+q2;
    %--- calculating dynamics of the system
    a1 = I1+M2*L1^2;
    a2 = I2+M2*D2^2;
    b = M2*L1*D2*cos(x(1)-x(2));
    D = a1*a2-b*b;

    xdot = zeros(4,1);
    xdot(1) = x(3); % angle in upper joint
    xdot(2) = x(4); % angle in lower joint
    xdot(3) = (a2*Q1-b*Q2)/D; % angular velocity in upper joint
    xdot(4) = (a1*Q2-b*Q1)/D; % angular velocity in lower joint
*/


//theta = [x(1), x(2)-x(1)];         %[shoulder elbow]
//thetaV = [x(3), x(4)-x(3)];        %[shoulder elbow]


typedef struct __biomech_data{
    //---- constructors
    __biomech_data( double m1, double l1, double m2, double l2 ) :
        M1( m1 ), L1( l1 ), M2( m2 ), L2( l2 ),
        BJ( 0.05 ), KJr( 3.0 ), BJr( 50.0 ),
        A1min( -pi ), A1max( pi ), A2min( -pi ), A2max( pi )
    {
        D1 = L1/2;
        D2 = L2/2;
        I1 = (M1*L1*L1)/3;
        I2 = (M2*L2*L2)/12;
        A1 = I1+M2*L1*L1;
        A2 = I2+M2*D2*D2;
        A12 = A1*A2;
        MLD = M2*L1*D2;
    };
    //---- methods
    void restriction( double a1min, double a1max, double a2min, double a2max )
    {
        //---- specify joint restrictions
        A1min = pi/2-a1max;
        A1max = pi/2-a1min;
        A2min = pi-a2max;
        A2max = pi-a2min;
    };

    //---- data
    // parameters of mechanical system
    double M1;    // mass of upper segment (kg)
    double M2;    // mass of lower segment (kg)
    double L1;    // lenght of upper segment (m)
    double L2;    // lenght of lower segment (m)

    double BJ;    // joint's viscosity
    double KJr;   // elasticity if joint's angle gets into restriction interval
    double BJr;   // viscosity if joint's angle gets into restriction interval

    // restriction angles for upper and lower joints
    double A1min; // minimal upper joint
    double A1max; // maximal upper joint
    double A2min; // minimal lower joint
    double A2max; // maximal lower joint

    // intermediated variables which are calculated in advance
    double D1; // L1/2; center mass for upper segment
    double D2; // L2/2; center mass for lower segment
    double I1; // (M1*L1*L1)/3; moment of inertia for upper segment
    double I2; // (M2*L2*L2)/12; moment of inertia for lower segment

    double A1;   // I1+M2*L1*L1
    double A2;   // I2+M2*D2*D2
    double A12;  // A1*A2
    double MLD;  // M2*L1*D2
} bm_data;

//---- calculating total non-conservative forces (F1n, F2n)
void calc_ncforces( void )
{
    //---- joints' angles and angular velocities
    double theta1 = x(1);       // shoulder
    double theta2 = x(2)-x(1);  // elbow
    double thetaV1 = x(3);      // shoulder
    double thetaV2 = x(4)-x(3); // elbow

    //--- joints viscosity
    double q1v = -BJ*thetaV(1);
    double q2v = -BJ*thetaV(2);
    //--- joints restrictions (elastic and viscosity components)
    double q1e = 0;
    double q2e = 0;
    if( theta1 > A1max && thetaV1 > 0 ){      // if shoulder flexion exceeds the limit
        q1e = -KJr*(theta1-A1max)-BJr*thetaV1;
    }
    else if( theta1 < A1min && thetaV1 < 0 ){ // if shoulder extension exceeds the limit
        q1e = -KJr*(theta1-A1min)-BJr*thetaV1;
    }
    if( theta2 > A2max && thetaV2 > 0 ){      // if elbow flexion exceeds the limit
        q2e = -KJr*(theta2-A2max)-BJr*thetaV2;
    }
    else if( theta2 < A2min && thetaV2 < 0 ){ // if elbow flexion exceeds the limit
        q2e = -KJr*(theta2-A2min)-BJr*thetaV2;
    }
    //---- calculating muscle force
/*
    [l, v, h_up, h_dw] = geometry_block( theta, thetaV ); %get muscle geomerty
    f = muscles( mn, l, v );
    double q1m = sum(f.*h_up);
    double q2m = sum(f.*h_dw);
*/
    //---- calculating total non-conservative forces
    F1n = q1v-q2v+q1e-q2e+q1m-q2m; // upper joint
    F2n = q2v+q2e+q2m;             // lower joint
}

//---- calculating total conservative forces (F1c, F2c)
void calc_cforces( void )
{
    double f1c = -MLD*sin(x(1)-x(2))*x(4)*x(4);
    double f2c = MLD*sin(x(1)-x(2))*x(3)*x(3);

    double D1g = D2g = 0;
    double f1g = M1*g*D1g*sin(x(1))+M2*g*L1g*sin(x(1));
    double f2g = M2*g*D2g*sin(x(2));

    F1c = f1c-f1g; // upper joint
    F2c = f2c-f2g; // lower joint
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
void solve_biomech( void )
{
    //---- calculating total conservative forces
    calc_cforces();  // get F1c and F2c forces
    //---- calculating total non-conservative forces
    calc_ncforces(); // get F1n and F2n forces

    double b = MLD*cos(x(1)-x(2));
    double d = A12-b*b;

    double q1 = F1c+F1n;      // total torque in upper joint;
    xdot(1) = x(3);           // angle in upper joint
    xdot(3) = (A2*q1-b*q2)/d; // angular velocity in upper joint

    double q2 = F2c+F2n;      // total torque in lower joint;
    xdot(2) = x(4);           // angle in lower joint
    xdot(4) = (A1*q2-b*q1)/d; // angular velocity in lower joint

/*
    Euler method
    k1 = step*F(t,y(t));
    y = y+k1;

    Runge-Kutta method (2nd order)
    k1 = step*F(t,y(t));
    k2 = step*F(t+h/2,y(t)+k1/2);
    y = y+k2;

    Runge-Kutta method (4th order)
    k1 = step*F(t,y(t))
    k2 = step*F(t+h/2,y(t)+k1/2)
    k3 = step*F(t+h/2,y(t)+k2/2)
    k4 = step*F(t+h,y(t)+k3)
    y = y+(k1+2*k2+2*k3+k4)/6;
*/
}

#endif // __BIOMECH_SKETCH_H
