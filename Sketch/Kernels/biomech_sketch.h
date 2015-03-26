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
        BJ( 0.05 ), KJr( 3.0 ), BJr( 50.0 )
    {
        M[0] = m1;
        L[0] = l1;
        D[0] = l1/2.;
        I[0] = ( m1*l1*l1 )/3.;


        M[1] = m2;
        L[1] = l2;
        D[1] = l2/2.;
        I[1] = ( m2*l2*l2 )/12.;

        ThetaMin[0] = -pi; ThetaMax[0] = pi;
        ThetaMin[1] = -pi; ThetaMax[1] = pi;

        A[0] = I[1]+M[1]*D[1]*D[1];
        A[1] = I[0]+M[1]*L[0]*L[0];
        A[2] = A[0]*A[1];
        MLD = M[1]*L[0]*D[1];
        L1g = D1g = D2g = 0;
    };

    //---- methods
    void restriction( double a1min, double a1max, double a2min, double a2max )
    {
        //---- specify joint restrictions
        ThetaMin[0] = pi/2-a1max;
        ThetaMax[0] = pi/2-a1min;
        ThetaMin[1] = pi-a2max;
        ThetaMax[1] = pi-a2min;
    };

    //---- data
    // parameters of mechanical system
    double M[2];  // mass of upper/lower segments (kg)
    double L[2];  // lenght of upper/lower segment (m)

    double BJ;    // joint's viscosity
    double KJr;   // elasticity if joint's angle gets into restriction interval
    double BJr;   // viscosity if joint's angle gets into restriction interval

    // restriction angles for upper and lower joints
    double ThetaMin[2]; // minimal angle for upper joint
    double ThetaMax[2]; // maximal angle for upper joint

    // intermediate variables which are calculated one time in advance of simulation
    double D[2];  // L1/2; center mass for upper/lower segment
    double I[2];  // (M1*L1*L1)/3; moment of inertia for upper segment; (M2*L2*L2)/12; moment of inertia for lower segment

    double L1g;   // projection L1 onto gravity vector
    double D1g;   // projection D1 onto gravity vector
    double D2g;   // projection D2 onto gravity vector

    double A[3];  // I[1]+M[1]*D[1]*D[1]; I[0]+M[1]*L[0]*L[0]; A[0]*A[1]
    double MLD;   // M[1]*L[0]*D[1]
} bm_data;

//---- calculating total non-conservative forces (Q1n, Q2n)
void calc_ncforces( double &q1n, double &q2n )
{
    //--- viscosity in upper joint
    double q1v = calc_vforce( thetaV[0] );
    //--- joints restrictions (elastic and viscosity components) for upper joint
    double q1e = calc_eforce( theta[0]-ThetaMax[0], thetaV[0] )*calc_eforce( -( theta[0]-ThetaMin[0] ), -thetaV[0] );

    //--- viscosity in lower joint
    double q2v = calc_vforce( thetaV[1] );
    //--- joints restrictions (elastic and viscosity components) for lower joint
    double q2e = calc_eforce( theta[1]-ThetaMax[1], thetaV[1] )*calc_eforce( -( theta[1]-ThetaMin[1] ), -thetaV[1] );

    //---- calculating muscle force
    double q1m = 0;
    double q2m = 0;
/*
    [l, v, h_up, h_dw] = geometry_block( theta, thetaV ); %get muscle geomerty
    f = muscles( mn, l, v );
    double q1m = sum(f.*h_up);
    double q2m = sum(f.*h_dw);
*/
    //---- calculating total non-conservative forces
    q1n = q1v+q1e+q1m-q2v-q2e-q2m; // upper joint
    q2n = q2v+q2e+q2m;             // lower joint
}

double calc_ncforce( double x, double v, double xmin, double xmax )
{
    double qv = calc_vforce( v );
    double qe = calc_eforce( x-xmax, v )*calc_eforce( -( x-xmin ), -v );
    return qv+qe;
}

double calc_vforce( double v )
{
    return -BJ*v;
}

double calc_eforce( double x, double v )
{
    int mul = int( x > 0 && v > 0 );
    return -( KJr*x+BJr*v )*mul;
}

//---- calculating total conservative forces (q1c, q2c)
double calc_cforce( double x1, double x2 )
{
    return MLD*sin(x1)*x2*x2;
/*
    double f1c = -MLD*sin(x(1)-x(2))*x(4)*x(4);
    double f2c = MLD*sin(x(1)-x(2))*x(3)*x(3);

    double f1g = -( M[0]*g*D1g*sin(x(1))+M[1]*g*L1g*sin(x(1)));
    double f2g = -( M[1]*g*D2g*sin(x(2)));

    q1c = f1c+f1g; // upper joint
    q2c = f2c+f2g; // lower joint
*/
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
    //---- joints' angles and angular velocities
    double theta[2] = {x(1s), x(2)-x(1)};  // shoulder; elbow
    double thetaV[2] = {x(3), x(4)-x(3)}; // shoulder; elbow
//  <synch> 1
    //---- pre-calculations
    double b = MLD*cos(x(1)-x(2));
    double d = A[2]-b*b;
    //---- calculating total conservative forces
    double q1c = calc_cforce( x(2)-x(1), x(4)); // upper segment
    double q2c = calc_cforce( x(1)-x(2), x(3)); // lower segment
    //---- calculating total non-conservative forces (excluding muscles)
    double q1n = calc_ncforce( theta[0], thetaV[0], ThetaMin[0], ThetaMax[0] ); // upper segment
    double q2n = calc_ncforce( theta[1], thetaV[1], ThetaMin[1], ThetaMax[1] ); // lower segment
    //---- calculating non-conservative forces produced by muscles
    double q1m = 0;
    double q2m = 0;
// <synch>
    //---- calculating total conservative and non-conservative forces
    double q1 = q1c+q1n+q1m-q2n-q2m; // total torque in upper joint;
    double q2 = q2c+q2n+q2m;          // total torque in lower joint;
// <synch>
    //---- solve the system of differential equations
    xdot(3) = (A[0]*q1-b*q2)/d;       // angular velocity in upper joint
    xdot(4) = (A[1]*q2-b*q1)/d;       // angular velocity in lower joint
// <synch>
    xdot(1) = x(3);                   // angle in upper joint
    xdot(2) = x(4);                   // angle in lower joint
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
