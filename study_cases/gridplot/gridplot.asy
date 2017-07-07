settings.outformat="pdf";
settings.prc=false;
settings.render=0;
import three;
currentprojection=perspective(0.5,0.3,0.1);


// Colors
pen spinup = rgb("226AB2");
pen spindown = rgb("B22222");
pen nicegray = rgb("AAAAAA");

real spinsize = 0.1;

size(7cm,7cm);

int L = 3;

rand();

srand(42);

for(int i=0; i<L; ++i){
  for(int j=0; j<L; ++j){
    for(int k=0; k<L; ++k){
      if(i>0) draw((i*X + j*Y + k*Z)--((i-1)*X + j*Y + k*Z),nicegray);
      if(j>0) draw((i*X + j*Y + k*Z)--(i*X + (j-1)*Y + k*Z),nicegray);
      if(k>0) draw((i*X + j*Y + k*Z)--(i*X + j*Y + (k-1)*Z),nicegray);
    }
  }
 }


for(int i=0; i<L; ++i){
  for(int j=0; j<L; ++j){
    for(int k=0; k<L; ++k){
      if (unitrand()> 0.5)
        draw((-spinsize*Z + i*X + j*Y + k*Z)--(spinsize*Z+ i*X + j*Y + k*Z)  , spinup, arrow=Arrow3());
      else
        draw((spinsize*Z + i*X + j*Y + k*Z)--(-spinsize*Z+ i*X + j*Y + k*Z)  , spindown, arrow=Arrow3());
    }
  }
 }

draw( (O + (L-1 + 0.3)*X) -- (O + (L-1 + 0.3)*X + (L-1)*Y)  , black,
      L=Label("L"
              , position=MidPoint
              , filltype=Fill(white)));
