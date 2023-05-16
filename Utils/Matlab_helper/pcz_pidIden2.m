function [W,C] = pcz_pidIden2(U,Y,Ts,args)
arguments
    U(:,1)
    Y(:,1)
    Ts(1,1)
    args.Type = 'PIDF'
end

    % Position output
    y    = Y(1:end-5);
    yp   = Y(2:end-4);
    ypp  = Y(3:end-3);
    yppp = Y(4:end-2);
    
    % Acceleration input
    u    = U(1:end-5);
    up   = U(2:end-4);
    upp  = U(3:end-3);
    uppp = U(4:end-2);
    
    % 2.-od rendu
    Reg_phi = [-yp -y upp up u];
    Reg_y = ypp;
    X = Reg_phi\Reg_y;
    az = [1 X(1:2)'];
    bz = X(3:5)';
    
    % 3.-ad rendu
    % Reg_phi = [-ypp -yp -y uppp upp up u];
    % Reg_y = yppp;
    % X = Reg_phi\Reg_y;
    % az = [1 X(1:3)'];
    % bz = X(4:7)';
    
    W = tf(bz,az,Ts);

    %{
    
    % Model order
    n = numel(az)-1;

    % Simulate (IIR filter)
    Ta = toeplitz(az(1:n),[1 zeros(n-1,1)]);
    Tb = toeplitz(bz(1:n),[bz(1) zeros(n-1,1)]);
    z0 = Ta*y(1:n) - Tb*u(1:n);
    y_filt = filter(bz,az,u,z0);

    % Simulate (state-space)
    S = ss(W);
    x0 = obsv(S) \ (y(1:n) - [S.D 0 ; S.C*S.B S.D]*u(1:2));
    y_lsim = lsim(S,u,[],x0);
        
    %}
    
    C = pidtune(W,args.Type);

end