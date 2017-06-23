settings.outformat="pdf";
settings.prc=false;
settings.render=0;
import three;
currentprojection=perspective(1,0.5,0.3);


// Colors
pen niceblue = rgb("3188A2");
pen niceyellow = rgb("8E7D42");
pen nicegray = rgb("767886");

real sep = 0.85;
real spinsize = 0.1;
real Jspan = 0.1;

size(5cm,5cm);

// Axis, let the J labels breathe.
// x
draw(-1*X--(-sep/2*X - Jspan*X),nicegray);
draw((-sep/2*X + Jspan*X)--O, nicegray);
draw(O--(sep/2*X - Jspan*X), nicegray);
draw((sep/2*X + Jspan*X)--X, nicegray);
// y
draw(-1*Y--(-sep/2*Y - Jspan*Y),nicegray);
draw((-sep/2*Y + Jspan*Y)--O, nicegray);
draw(O--(sep/2*Y - Jspan*Y), nicegray);
draw((sep/2*Y + Jspan*Y)--Y, nicegray);
// z
draw(-1*Z--(-sep/2*Z - Jspan*Z),nicegray);
draw((-sep/2*Z + Jspan*Z)--O, nicegray);
draw(O--(sep/2*Z - Jspan*Z), nicegray);
draw((sep/2*Z + Jspan*Z)--Z, nicegray);


// Spins
draw((-spinsize*Z)--(spinsize*Z)  , niceblue, arrow=Arrow3());

draw((-spinsize*Z - sep*X)--(spinsize*Z - sep*X)  , niceblue, arrow=Arrow3());
draw((-spinsize*Z + sep*X)--(spinsize*Z + sep*X)  , niceblue, arrow=Arrow3());

draw((+spinsize*Z - sep*Y)--(-spinsize*Z - sep*Y) , niceblue, arrow=Arrow3());
draw((-spinsize*Z + sep*Y)--(spinsize*Z + sep*Y)  , niceblue, arrow=Arrow3());

draw((+spinsize*Z - sep*Z)--(-spinsize*Z - sep*Z) , niceblue, arrow=Arrow3());
draw((+spinsize*Z + sep*Z)--(-spinsize*Z + sep*Z) , niceblue, arrow=Arrow3());

// J's
label(Label("+", niceyellow, filltype=Fill(white)), +sep/2*X);
label(Label("+", niceyellow, filltype=Fill(white)), -sep/2*X);
label(Label("-", niceyellow, filltype=Fill(white)), +sep/2*Y);
label(Label("-", niceyellow, filltype=Fill(white)), -sep/2*Y);
label(Label("+", niceyellow, filltype=Fill(white)), +sep/2*Z);
label(Label("+", niceyellow, filltype=Fill(white)), -sep/2*Z);
