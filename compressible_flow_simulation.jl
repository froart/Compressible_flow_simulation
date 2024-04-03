using GLMakie, GeometryBasics, Observables, LinearAlgebra

GLMakie.activate!()

fig = Figure(size = (600, 600))
ax  = Axis3(fig[1,1][1,1], width = 500, height = 500)
# Quit upon ESC key
quit = Observable(false)
on(events(fig).keyboardbutton) do event
   if event.action == Keyboard.press && event.key == Keyboard.escape
      quit[] = true
      notify(quit)
   end
end

# defining initial valUes 
N2 = 31;
M2 = 31;

U = zeros(N2, M2, 2); # velocities
E = zeros(N2, M2); # energy
RO = zeros(N2, M2); # density
P = zeros(N2, M2); # pressure
US = zeros(N2, M2, 2); # velocity
ES = zeros(N2, M2); # ?
ROS = zeros(N2, M2); # ?
MW = zeros(N2, M2); # ?

Max = 2.0; # Mach nUmber
g = 1.4; # adiabata nUmber
dx = 0.1; # grid x step
dy = dx; # grid y step
t = 0.0; # cUrrent time 
tk = 10.0; # period of observation, ini: 20
dt = 0.1*dx; # step of the time
x=0:dx:(N2-1)*dx;
y=0:dy:(M2-1)*dy;

#### MOD1 ####
# defining intial field. HERE PARAMETERS HAVE TO CHANGE TO CREATE AN EXPLOSION.
for J in 1:M2
    for I in 1:N2
      RO[I,J]=1.0;
      U[I,J,1]=0.0;
      U[I,J,2]=0.0;
      P[I,J]=1.0;
      E[I,J]=1.0/((g-1)*1.0);
    end
end

# plot initial values
surface!(ax, x, y, RO);
zlims!(ax, (0.5, 1.5)) 


# Animation function
function step!()
    #### MOD2 ####
    for J = 1:M2
      RO[1,J]=RO[2,J];
      U[1,J,1]=U[2,J,1];
      U[1,J,2]=U[2,J,2];
      P[1,J]=P[2,J];
      E[1,J]=E[2,J];
      # правая граница:
      U[N2,J,1]=U[N2-1,J,1];
      U[N2,J,2]=U[N2-1,J,2];
      RO[N2,J]=RO[N2-1,J];
      E[N2,J]=E[N2-1,J];
    end
    # нижняя граница:
    for I=1:N2
      U[I,1,1]=U[I,2,1];
      U[I,1,2]=-U[I,2,2];
      RO[I,1]=RO[I,2];
      E[I,1]=E[I,2];
      # верхняя граница:
      U[I,M2,1]=U[I,M2-1,1];
      U[I,M2,2]=U[I,M2-1,2];
      RO[I,M2]=RO[I,M2-1];
      E[I,M2]=E[I,M2-1];
    end
 
  #### MOD4 ####
  # compUting the pressUre
    for J in 2:M2-1
      for I in 2:N2-1
        P[I,J]=RO[I,J]*(g-1).*(E[I,J]-.5*(U[I,J,1].^2+U[I,J,2].^2));
      end
    end
  #### MOD5 ####
    # compUting Euleur's stage
    for J in 2:M2-1
      for I in 2:N2-1
        P0=P[I,J];
        PL=.5*(P0+P[I-1,J]);
        PP=.5*(P0+P[I+1,J]);    
        PN=.5*(P0+P[I,J-1]);
        PV=.5*(P0+P[I,J+1]);
        U01=U[I,J,1];
        U02=U[I,J,2];
        UL1=.5*(U01+U[I-1,J,1]);
        UP1=.5*(U01+U[I+1,J,1]);
        UN2=.5*(U02+U[I,J-1,2]);
        UV2=.5*(U02+U[I,J+1,2]);
        # computing intermediate values of velocities and total energy
        US[I,J,1]=U[I,J,1]+(PL-PP)*dt/(dx*RO[I,J]);
        US[I,J,2]=U[I,J,2]+(PN-PV)*dt/(dy*RO[I,J]);
        ES[I,J]=E[I,J]+(PL*UL1-PP*UP1)*dt/(dx*RO[I,J])+(PN*UN2-PV*UV2)*dt/(dy*RO[I,J]);
      end
    end
  #### MOD2S ####
  # defining border conditions before Lagrange stage
    for J = 1:M2
      # левая граница: 
      ROS[1,J]=ROS[2,J];
      US[1,J,1]=US[2,J,1];
      US[1,J,2]=US[2,J,2];
      ES[1,J]=ES[2,J];
      # правая граница: 
      US[N2,J,1]=US[N2-1,J,1];
      US[N2,J,2]=US[N2-1,J,2];
      ES[N2,J]=ES[N2-1,J];
    end
    for I=1:N2
      # нижняя граница:
      US[I,1,1]=US[I,2,1];
      US[I,1,2]=-US[I,2,2];
      ES[I,1]=ES[I,2];
      # верхняя граница:
      US[I,M2,1]=US[I,M2-1,1];
      US[I,M2,2]=US[I,M2-1,2];
      ES[I,M2]=ES[I,M2-1];
    end

  #### MOD10 ####
  # Lagrange and final stages
    for I in 2:N2-1;
      for J in 2:M2-1;
        RO0=RO[I,J];
        U10=US[I,J,1];
        U20=US[I,J,2];
        UL=.5*(U10+US[I-1,J,1]);
        UP=.5*(U10+US[I+1,J,1]);
        UN=.5*(U20+US[I,J-1,2]);
        UV=.5*(U20+US[I,J+1,2]);
        # defining direction of the velocity and computing flows 
        # of mass throught the left border of the cell
        if UL <= 0.0
          DM1=UL*RO0*dy*dt;
          D1=0.0;
        else 
          DM1=UL*RO[I-1,J]*dy*dt;
          D1=1.0;
        end
        # defining direction of the velocity and computing flows 
        # of mass throught the lower border of the cell
        if UN <= 0.0
          DM2=UN*RO0*dx*dt;
          D2=0.0;
        else 
          DM2=UN*RO[I,J-1]*dx*dt;
          D2=1.0;
        end
        # defining direction of the velocity and computing flows 
        # of mass throught the right border of the cell
        if UP <= 0.0
          DM3=UP*RO[I+1,J]*dy*dt;
          D3=1.0;
        else
          DM3=UP*RO0*dy*dt;
          D3=0.0;
        end
        # defining direction of the velocity and computing flows 
        # of mass throught the upper border of the cell
        if UV <= 0.0
          DM4=UV*RO[I,J+1]*dx*dt;
          D4=1.0;
        else
          DM4=UV*RO0*dx*dt;
          D4=0.0;
        end
        
        AM1=abs(DM1);
        AM2=abs(DM2);
        AM3=abs(DM3);
        AM4=abs(DM4);
        Z1=dy*dx;
        Z2=(D1-0.5)*AM1+(D2-0.5)*AM2+(D3-0.5)*AM3+(D4-0.5)*AM4;
        # computing the final value of the density
        ROS[I,J]=RO0+2.0*Z2/Z1;
        Z3=RO0*Z1-(1.0-D1)*AM1-(1.0-D2)*AM2-(1.0-D3)*AM3-(1.0-D4)*AM4;
        Z4=ROS[I,J]*Z1;
        B1=D1*AM1;
        B2=D2*AM2;
        B3=D3*AM3;
        B4=D4*AM4;
        # computing final values of the components of the speed
        U[I,J,1]=(U10*Z3+US[I-1,J,1]*B1+US[I,J-1,1]*B2+US[I+1,J,1]*B3+US[I,J+1,1]*B4)/Z4;
        U[I,J,2]=(U20*Z3+US[I-1,J,2]*B1+US[I,J-1,2]*B2+US[I+1,J,2]*B3+US[I,J+1,2]*B4)/Z4;
        # computing the final value of the total energy
        E[I,J]=(ES[I,J]*Z3+ES[I-1,J]*B1+ES[I,J-1]*B2+ES[I+1,J]*B3+ES[I,J+1]*B4)/Z4;
      end
    end

  # restore the true values of the density
  for I in 2:N2-1
    for J in 2:M2-1
      RO[I,J]=ROS[I,J];
    end
  end

  return x, y, RO
end 

# Run the animation
fps = 60
frame_count = 0

# Boom here
RO[floor(Int, N2/2), floor(Int, M2/2)] = 15.0

while !quit[]
    display(fig)
    local x, y, RO = step!()
    empty!(ax)
    surface!(ax, x, y, RO);
    sleep(1/fps)
    global frame_count = frame_count + 1
end

GLMakie.closeall()
