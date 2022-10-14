function [z,rate,height,width]=GaragashAndGermanovich2022Arxiv(Kc,E,nu,eta,deltagamma,t,V)


    %Garagash and Germanovich Aug 2022 Arxiv - Page 16
    Kprime=sqrt(2/pi)*Kc;
    Eprime=(1/pi)*(E/(1-nu^2));
    muprime=pi^2*eta;
    Vg=0.408*(Kprime^(8/3)/(Eprime*(deltagamma)^(5/3)));
    z2=2.2*(((V-Vg).^2.*deltagamma.^(7/3).*t)./(muprime*Kprime^(4/3))).^(1/3);
    rate2=gradient(z2)./gradient(t);
    z=z2;rate=rate2;

    height=(Kprime/(deltagamma))^(2/3);
    width=height*0.396;

end