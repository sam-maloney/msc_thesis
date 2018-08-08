pOwn = 2.26853428554387;
UOwn = [0 -0.172084298482207 0];
TOwn = 1.38471475550203;
ROwn = 0.714285714285714;
CvOwn = 2.14285714285714;
phiOwn = 0.0427282724469507;

pNei = 2.23325954860948;
UNei = [0 -0.174118876651778 0];
TNei = 1.37930014895722;
RNei = 0.714285714285714;
CvNei = 2.14285714285714;
phiNei = 0.0578463881495107;

pLeft = 2.27525153216717;
ULeft = [0 -0.171703123115643 0];
TLeft = 1.38573867214271;

phiFace = 0.0502873302982307;

pRight = 2.25543910857366;
URight = [0 -0.172833102608523 0];
TRight = 1.38271208903968;

gammaOwn = ROwn/CvOwn+1;
rhoOwn = pOwn/(ROwn*TOwn);
KOwn = pOwn/rhoOwn^gammaOwn;
mOwn = rhoOwn*UOwn;
m2Own = norm(mOwn)^2;
rhoLeft = rhoOwn;
rhoPrev = -1;

while abs(rhoLeft - rhoPrev) > 2e-15
    rhoPrev = rhoLeft;
    rhoLeft = rhoLeft...
        - (0.5*m2Own/rhoLeft^2 + gammaOwn/(gammaOwn-1)*KOwn*rhoLeft^(gammaOwn-1) + phiFace - B)...
        / (-m2Own/rhoLeft^3 + gammaOwn*KOwn*rhoLeft^(gammaOwn-2));
end