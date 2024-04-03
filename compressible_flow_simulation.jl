using Base, GLMakie, GeometryBasics, Observables, LinearAlgebra

Base.@kwdef mutable struct LargeParticleMethod
    N::Int64 = 30
    M::Int64 = 30
    v::Array{Float64, 3} # velocity 
    E::Array{Float64, 2} # energy
    ρ::Array{Float64, 2} # density
    P::Array{Float64, 2} # pressure
    v_n::Array{Float64, 3} # new value of velocity 
    E_n::Array{Float64, 2} # new value of energy
    ρ_n::Array{Float64, 2} # new value of density
    dx::Float64 = 0.1
    dy::Float64 = 0.1
    dt::Float64 = 0.01
    x
    y
    g::Float64  = 1.4 # adiabata number
end

# Animation function
function step!(lpm::LargeParticleMethod)
    #### MOD2 ####
    for j in 1:lpm.M
        lpm.ρ[1,j] = lpm.ρ[2,j];
        lpm.v[1,j,1]=lpm.v[2,j,1];
        lpm.v[1,j,2]=lpm.v[2,j,2];
        lpm.P[1,j]=lpm.P[2,j];
        lpm.E[1,j]=lpm.E[2,j];
        # правая граница:
        lpm.v[lpm.N,j,1]=lpm.v[lpm.N-1,j,1];
        lpm.v[lpm.N,j,2]=lpm.v[lpm.N-1,j,2];
        lpm.ρ[lpm.N,j]=lpm.ρ[lpm.N-1,j];
        lpm.E[lpm.N,j]=lpm.E[lpm.N-1,j];
    end
    # нижняя граница:
    for i=1:lpm.N
      lpm.v[i,1,1]=lpm.v[i,2,1];
      lpm.v[i,1,2]=-lpm.v[i,2,2];
      lpm.ρ[i,1]=lpm.ρ[i,2];
      lpm.E[i,1]=lpm.E[i,2];
      # верхняя граница:
      lpm.v[i,lpm.M,1]=lpm.v[i,lpm.M-1,1];
      lpm.v[i,lpm.M,2]=lpm.v[i,lpm.M-1,2];
      lpm.ρ[i,lpm.M]=lpm.ρ[i,lpm.M-1];
      lpm.E[i,lpm.M]=lpm.E[i,lpm.M-1];
    end
 
  #### MOD4 ####
  # computing the pressure
    for j in 2:lpm.M-1
      for i in 2:lpm.N-1
        lpm.P[i,j]=lpm.ρ[i,j]*(lpm.g-1).*(lpm.E[i,j]-.5*(lpm.v[i,j,1].^2+lpm.v[i,j,2].^2));
      end
    end
  #### MOD5 ####
    # computing Euler's stage
    for j in 2:lpm.M-1
      for i in 2:lpm.N-1
        P0=lpm.P[i,j];
        PL=.5*(P0+lpm.P[i-1,j]);
        PP=.5*(P0+lpm.P[i+1,j]);    
        PN=.5*(P0+lpm.P[i,j-1]);
        PV=.5*(P0+lpm.P[i,j+1]);
        U01=lpm.v[i,j,1];
        U02=lpm.v[i,j,2];
        UL1=.5*(U01+lpm.v[i-1,j,1]);
        UP1=.5*(U01+lpm.v[i+1,j,1]);
        UN=.5*(U02+lpm.v[i,j-1,2]);
        UV2=.5*(U02+lpm.v[i,j+1,2]);
        # computing intermediate values of velocities and total energy
        lpm.v_n[i,j,1]=lpm.v[i,j,1]+(PL-PP)*lpm.dt/(lpm.dx*lpm.ρ[i,j]);
        lpm.v_n[i,j,2]=lpm.v[i,j,2]+(PN-PV)*lpm.dt/(lpm.dy*lpm.ρ[i,j]);
        lpm.E_n[i,j]=lpm.E[i,j]+(PL*UL1-PP*UP1)*lpm.dt/(lpm.dx*lpm.ρ[i,j])+(PN*UN-PV*UV2)*lpm.dt/(lpm.dy*lpm.ρ[i,j]);
      end
    end
  #### MOD2S ####
  # defining border conditions before Lagrange stage
    for j = 1:lpm.M
      # левая граница: 
      lpm.ρ[1,j]=lpm.ρ[2,j];
      lpm.v_n[1,j,1]=lpm.v_n[2,j,1];
      lpm.v_n[1,j,2]=lpm.v_n[2,j,2];
      lpm.E_n[1,j]=lpm.E_n[2,j];
      # правая граница: 
      lpm.v_n[lpm.N,j,1]=lpm.v_n[lpm.N-1,j,1];
      lpm.v_n[lpm.N,j,2]=lpm.v_n[lpm.N-1,j,2];
      lpm.E_n[lpm.N,j]=lpm.E_n[lpm.N-1,j];
    end
    for i=1:lpm.N
      # нижняя граница:
      lpm.v_n[i,1,1]=lpm.v_n[i,2,1];
      lpm.v_n[i,1,2]=-lpm.v_n[i,2,2];
      lpm.E_n[i,1]=lpm.E_n[i,2];
      # верхняя граница:
      lpm.v_n[i,lpm.M,1]=lpm.v_n[i,lpm.M-1,1];
      lpm.v_n[i,lpm.M,2]=lpm.v_n[i,lpm.M-1,2];
      lpm.E_n[i,lpm.M]=lpm.E_n[i,lpm.M-1];
    end

  #### MOD10 ####
  # Lagrange and final stages
    for i in 2:lpm.N-1;
      for j in 2:lpm.M-1;
        ρ0=lpm.ρ[i,j];
        U10=lpm.v_n[i,j,1];
        U20=lpm.v_n[i,j,2];
        UL=.5*(U10+lpm.v_n[i-1,j,1]);
        UP=.5*(U10+lpm.v_n[i+1,j,1]);
        UN=.5*(U20+lpm.v_n[i,j-1,2]);
        UV=.5*(U20+lpm.v_n[i,j+1,2]);
        # defining direction of the velocity and computing flows 
        # of mass through the left border of the cell
        if UL <= 0.0
          DM1=UL*ρ0*lpm.dy*lpm.dt;
          D1=0.0;
        else 
          DM1=UL*lpm.ρ[i-1,j]*lpm.dy*lpm.dt;
          D1=1.0;
        end
        # defining direction of the velocity and computing flows 
        # of mass through the lower border of the cell
        if UN <= 0.0
          DM=UN*ρ0*lpm.dx*lpm.dt;
          D2=0.0;
        else 
          DM=UN*lpm.ρ[i,j-1]*lpm.dx*lpm.dt;
          D2=1.0;
        end
        # defining direction of the velocity and computing flows 
        # of mass through the right border of the cell
        if UP <= 0.0
          DM3=UP*lpm.ρ[i+1,j]*lpm.dy*lpm.dt;
          D3=1.0;
        else
          DM3=UP*ρ0*lpm.dy*lpm.dt;
          D3=0.0;
        end
        # defining direction of the velocity and computing flows 
        # of mass through the upper border of the cell
        if UV <= 0.0
          DM4=UV*lpm.ρ[i,j+1]*lpm.dx*lpm.dt;
          D4=1.0;
        else
          DM4=UV*ρ0*lpm.dx*lpm.dt;
          D4=0.0;
        end
        
        AM1=abs(DM1);
        AM=abs(DM);
        AM3=abs(DM3);
        AM4=abs(DM4);
        Z1=lpm.dy*lpm.dx;
        Z2=(D1-0.5)*AM1+(D2-0.5)*AM+(D3-0.5)*AM3+(D4-0.5)*AM4;
        # computing the final value of the density
        lpm.ρ_n[i,j]=ρ0+2.0*Z2/Z1;
        Z3=ρ0*Z1-(1.0-D1)*AM1-(1.0-D2)*AM-(1.0-D3)*AM3-(1.0-D4)*AM4;
        Z4=lpm.ρ_n[i,j]*Z1;
        B1=D1*AM1;
        B2=D2*AM;
        B3=D3*AM3;
        B4=D4*AM4;
        # computing final values of the components of the speed
        lpm.v[i,j,1]=(U10*Z3+lpm.v_n[i-1,j,1]*B1+lpm.v_n[i,j-1,1]*B2+lpm.v_n[i+1,j,1]*B3+lpm.v_n[i,j+1,1]*B4)/Z4;
        lpm.v[i,j,2]=(U20*Z3+lpm.v_n[i-1,j,2]*B1+lpm.v_n[i,j-1,2]*B2+lpm.v_n[i+1,j,2]*B3+lpm.v_n[i,j+1,2]*B4)/Z4;
        # computing the final value of the total energy
        lpm.E[i,j]=(lpm.E_n[i,j]*Z3+lpm.E_n[i-1,j]*B1+lpm.E_n[i,j-1]*B2+lpm.E_n[i+1,j]*B3+lpm.E_n[i,j+1]*B4)/Z4;
      end
    end

  # restore the true values of the density
  for i in 2:lpm.N-1
    for j in 2:lpm.M-1
      lpm.ρ[i,j]=lpm.ρ_n[i,j];
    end
  end

  return lpm.x, lpm.y, lpm.ρ
end

function run_simulation(N, M, dx, dt)

    GLMakie.activate!()
    # Create plot
    fig = Figure(size = (700, 600))
    ax  = Axis3(fig[1,1][1,1], width = 500, height = 500, zlabel = "ρ (density)", xlabel = "x", ylabel = "y")
    zlims!(ax, (0.5, 1.5)) 
    # Quit upon ESC key
    quit = Observable(false)
    on(events(fig).keyboardbutton) do event
       if event.action == Keyboard.press && event.key == Keyboard.escape
          quit[] = true
          notify(quit)
       end
    end

    lpm = LargeParticleMethod(N   = N, 
                              M   = M,
                              v   = zeros(N, M, 2), 
                              E   = zeros(N, M),
                              ρ   = zeros(N, M),
                              P   = zeros(N, M),
                              v_n = zeros(N, M, 2),
                              E_n = zeros(N, M),
                              ρ_n = zeros(N, M),
                              dx  = dx,
                              dy  = dx,
                              x   = 0:dx:(N-1)*dx,
                              y   = 0:dx:(M-1)*dx,
                              dt  = dt) 

    # Intializing fields
    for j in 1:M, i in 1:N
        lpm.ρ[i,j]   = 1.0;
        lpm.v[i,j,1] = 0.0;
        lpm.v[i,j,2] = 0.0;
        lpm.P[i,j]   = 1.0;
        lpm.E[i,j]   = 1.0/((lpm.g-1)*1.0);
    end

    # Boom here
    x_c, y_c = floor(Int, N/2), floor(Int, M/2)
    lpm.ρ[x_c, y_c] = 15.0

    # Run the simulation
    fps = 60
    while !quit[]
        display(fig)
        local x, y, ρ = step!(lpm)
        empty!(ax)
        surface!(ax, x, y, ρ);
        sleep(1/fps)
    end

    GLMakie.closeall()
end
