%Topology File
tFileID = fopen('ieee57.cdf');
line = fgetl(tFileID);
line = fgetl(tFileID);
nBus = strsplit(line);
nBus = str2double(nBus{4});
fgetl(tFileID);
count999 = 0;
yBus = zeros(nBus, nBus);
yShunt = zeros(nBus, nBus);
V = ones(nBus,1); % Initialize the bus voltages..
theta = zeros(nBus,1); % Initialize the bus angles..
E = [theta(1:end); V];   % State Vector..
while count999<2
   line = fgetl(tFileID);
   if line(1:4) == '-999'
       count999 = count999 + 1;
       line = fgetl(tFileID);
       line = fgetl(tFileID);
   end

   if count999 == 0
        bus = str2double(line(1:4));
        yBus(bus, bus) = str2double(line(107:114)) + 1i*str2double(line(115:122));

   elseif count999 == 1
        i = str2double(line(1:4));
        j = str2double(line(6:9));

        resistance          = str2double(line(20:29));
        inductance          = str2double(line(30:40));
        chargingSusceptance = 1i*str2double(line(41:50));
        yseries             = 1/(resistance + 1i*inductance);

        a = str2double(line(77:82));
        if a == 0
            a = 1;
        end
        
        yBus(i, i) = yseries/a^2 + chargingSusceptance/2;
        yBus(i, j) = -yseries/a;
        yBus(j, i) = -yseries/a;
        yBus(j, j) = yseries + chargingSusceptance/2;
            
        yShunt(i, j) = yShunt(i, j) + chargingSusceptance/2 + yseries*(1-a)/a^2;
        yShunt(j, i) = yShunt(j, i) +chargingSusceptance/2 + yseries*(a-1)/a;
   end
end
G = real(yBus);
B = imag(yBus);