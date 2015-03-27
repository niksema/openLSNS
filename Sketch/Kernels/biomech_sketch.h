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

function [L, V, Hup, Hdw] = geometry_block( theta, thetaV )
global A1min A1max A2min A2max;    % restriction angles for upper and lower joints
global Lopt R

% muscles map
% [shoulder_flexor shoulder_extensor elbow_flexor elbow_extensor
% bifunc_flexor bifunc_extensor]

    lm = 0.97*[A1max-A1min, A2max-A2min];                 %range
    phi = [pi/2-theta(1), pi-theta(2)];
    alpha = [phi(1)-A1min, A1max-phi(1), phi(2)-A2min, A2max-phi(2)];
    alpha = max( alpha, 0 );
    a = [alpha(1)/lm(1), alpha(2)/lm(1), alpha(3)/lm(2), alpha(4)/lm(2), ...
        (alpha(1)+alpha(3))/(lm(1)+lm(2)), (alpha(2)+alpha(4))/(lm(1)+lm(2))];
    L = a.*Lopt;                                      %muscle lenghts
    vA = [-thetaV(1), -thetaV(2)];                    %positive then extension, negative then flexion
    V = [vA(1)*R(1), vA(1)*R(2), vA(2)*R(3), vA(2)*R(4), ...
         vA(1)*R(5)+vA(2)*R(7), vA(1)*R(6)+vA(2)*R(8)]; %muscles velocities
    Hup = [R(1), R(2), 0,    0,    R(5), R(6)];
    Hdw = [0,    0,    R(3), R(4), R(7), R(8)];

% Calculate muscle forces
% MN - motoneuronal activation
% L - lenght of muscles
% V - muscles velocity
function F = muscles( MN, L, V )
    global Fmax Lopt Beta Omega Ro Av1 Bv1 Bv2; % properties of the muscle model
    lnorm = max( L./Lopt, 0 ); %normalized lenght
    Fl_active = exp( -(abs( (lnorm.^Beta-1))/Omega ).^Ro );
    FL_passive = 3.5*log( exp( (lnorm-1.4)./0.05 )+1.0 )...
                - 0.02*( exp( -18.7*(lnorm-0.79) )-1.0 );
    FL_passive = FL_passive.*max( lnorm-1, 0.0 );
    FV = zeros( 1, length(V) );
    for i = 1:length(V)
        if V(i) <= 0.0
            FV(i) = (Bv1-Av1*V(i))/(V(i)+Bv1);
        else
            FV(i) = (Bv2-(-5.34*lnorm(i)*lnorm(i)...
                   + 8.41* lnorm(i)-4.7)*V(i))/(V(i)+Bv2);
        end
    end;
    F = Fmax.*( MN.*Fl_active.*FV+FL_passive );


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

//---- calculating Coriolis conservative force
double calc_cforce( double x, double v )
{
    return MLD*sin(x)*v*v;
}

//---- calculating gravity conservative forces
double calc_gforce( double m, double l, double x )
{
    return m*g*l*sin(x);
}

//---- calculating non-conservative viscosity force
double calc_vforce( double v )
{
    return -BJ*v;
}

//---- calculating non-conservative elastic force
double calc_eforce( double x, double v )
{
    int mul = int( x > 0 && v > 0 );
    return -( KJr*x+BJr*v )*mul; // return ( x > 0 && v > 0 )? -KJr*x-BJr*v: 0;
}

//---- calculating non-conservative muscle force
double calc_mforce( double mn, double lnorm, double v )
{
    lnorm = (lnorm>0.)? lnorm:0.;
    double fl_a = exp(-(abs((lnorm^Beta-1))/Omega )^Ro ); // active component
    double fl_p = 3.5*log(exp((lnorm-1.4)/0.05 )+1.0 )-0.02*(exp(-18.7*(lnorm-0.79))-1.0); // passive component
    fl_p = (fl_p>0.)? fl_p:0.;
    fv = (v>0.)?(Bv2-(-5.34*lnorm*lnorm+8.41*lnorm-4.7)*v)/(v+Bv2):(Bv1-Av1*v)/(v+Bv1);    // velocity component
    return Fmax*( mn*fl_a*fv+fl_p);
}

/*
function [L, V, Hup, Hdw] = geometry_block( theta, thetaV )
global A1min A1max A2min A2max;    % restriction angles for upper and lower joints
global Lopt R

% muscles map
% [shoulder_flexor shoulder_extensor elbow_flexor elbow_extensor
% bifunc_flexor bifunc_extensor]

    lm = 0.97*[A1max-A1min, A2max-A2min];                 %range
    phi = [pi/2-theta(1), pi-theta(2)];
    alpha = [phi(1)-A1min, A1max-phi(1), phi(2)-A2min, A2max-phi(2)];
    alpha = max( alpha, 0 );
    a = [alpha(1)/lm(1), alpha(2)/lm(1), alpha(3)/lm(2), alpha(4)/lm(2), ...
        (alpha(1)+alpha(3))/(lm(1)+lm(2)), (alpha(2)+alpha(4))/(lm(1)+lm(2))];
    L = a.*Lopt;                                      %muscle lenghts
    Hup = [R(1), R(2), 0,    0,    R(5), R(6)];
    Hdw = [0,    0,    R(3), R(4), R(7), R(8)];

    vA = [-thetaV(1), -thetaV(2)];                    %positive then extension, negative then flexion
    V = [vA(1)*R(1), vA(1)*R(2), vA(2)*R(3), vA(2)*R(4), ...
         vA(1)*R(5)+vA(2)*R(7), vA(1)*R(6)+vA(2)*R(8)]; %muscles velocities
*/

double get_musclel( double x_up, double x_dw )
{
    return 0;
/*
    lm = 0.97*[A1max-A1min, A2max-A2min];                <-- calculate in advance

    phi = [pi/2-theta(1), pi-theta(2)];
    alpha = [phi(1)-A1min, A1max-phi(1), phi(2)-A2min, A2max-phi(2)];
    alpha = max( alpha, 0 );
    a = [alpha(1)/lm(1), alpha(2)/lm(1), alpha(3)/lm(2), alpha(4)/lm(2), ...
        (alpha(1)+alpha(3))/(lm(1)+lm(2)), (alpha(2)+alpha(4))/(lm(1)+lm(2))];
    L = a.*Lopt;                                      %muscle lenghts
    Hup = [R(1), R(2), 0,    0,    R(5), R(6)];
    Hdw = [0,    0,    R(3), R(4), R(7), R(8)];
*/
}

double get_musclev( double v_up, double v_dw )
{
    return 0;
/*
    vA = [-v_up, -v_dw];                    %positive then extension, negative then flexion
    V = [vA(1)*R(1), vA(1)*R(2), vA(2)*R(3), vA(2)*R(4), ...
         vA(1)*R(5)+vA(2)*R(7), vA(1)*R(6)+vA(2)*R(8)]; %muscles velocities
    Hup = [R(1), R(2), 0,    0,    R(5), R(6)];
    Hdw = [0,    0,    R(3), R(4), R(7), R(8)];
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
    double theta[2] = {x(1), x(2)-x(1)};         // shoulder; elbow
    double thetaV[2] = {x(3), x(4)-x(3)};        // shoulder; elbow
    /*
    for( int i = 0; i < 6; ++i ){
        lm[i] = get_musclel( theta[0], theta[1] );   // get muscle geometry
        vm[i] = get_musclev( thetaV[0], thetaV[1] ); // get muscle velocity
    }
    */
    double q1m[6] = {0};
    double q2m[6] = {0};
//  <synch> 1
    //---- pre-calculations
    double b = MLD*cos( x(1)-x(2) );
    double d = A[2]-b*b;
    //---- calculating Coriolis conservative forces
    double q1c = calc_cforce(x(2)-x(1),x(4));    // upper segment; double q1c = -MLD*sin(x(1)-x(2))*x(4)*x(4);
    double q2c = calc_cforce(x(1)-x(2),x(3));    // lower segment; double q2c = MLD*sin(x(1)-x(2))*x(3)*x(3);
    //---- calculating gravity conservative forces
    double q1g = -calc_gforce(M[0],D1g,x(1))-calc_gforce(M[1],L1g,x(1)); // upper segment; double q1g = -M[0]*g*D1g*sin(x(1))-M[1]*g*L1g*sin(x(1));
    double q2g = -calc_gforce(M[1],D2g,x(2));                            // lower segment; double q2g = -M[1]*g*D2g*sin(x(2));
    //---- calculating total non-conservative forces (excluding muscles)
    double q1v = calc_vforce(thetaV[0]);         // viscosity for upper segment
    double q2v = calc_vforce(thetaV[1]);         // viscosity for lower segment
    double q1e1 = calc_eforce(theta[0]-ThetaMax[0],thetaV[0]);     // max restriction for upper segment
    double q1e2 = calc_eforce(-(theta[0]-ThetaMin[0]),-thetaV[0]); // min restriction for upper segment
    double q2e1 = calc_eforce(theta[1]-ThetaMax[1],thetaV[1]);     // max restriction for lower segment
    double q2e2 = calc_eforce(-(theta[1]-ThetaMin[1]),-thetaV[1]); // min restriction for lower segment
    //---- calculating non-conservative forces produced by muscles
    for( int i = 0; i < 6; ++i ){
    //<one thread>
        f[i] = calc_mforce( mn[i], lm[i]/lopt[i], vm[i] );
        q1m[i] = f[i]*h_up[i];
        q2m[i] = f[i]*h_dw[i];
    //</one thread>
    }
// <synch> 2
    //---- calculating total conservative and non-conservative forces
    double q1 = q1c+q1g+q1v+q1e1+q1e2+q1m-q2v-q2e1-q2e2-q2m; // total torque in upper joint;
    double q2 = q2c+q2g+q2v+q2e1+q2e2+q2m;                   // total torque in lower joint;
// <synch> 3
    //---- solve the system of differential equations
    xdot(3) = (A[0]*q1-b*q2)/d;                  // angular velocity in upper joint
    xdot(4) = (A[1]*q2-b*q1)/d;                  // angular velocity in lower joint
    xdot(1) = x(3);                              // angle in upper joint
    xdot(2) = x(4);                              // angle in lower joint
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
