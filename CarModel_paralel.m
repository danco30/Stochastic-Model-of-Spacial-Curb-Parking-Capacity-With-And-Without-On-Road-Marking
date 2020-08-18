% curb parking
% three kinds of spaces:
% white - available
% gray - parking
% black - obstile
clear all
clc

format longg
rng(123); % creating a single seed for repetetivity

% Define street
L_street = 150; % [meters]
L_street_binning = L_street*100; % [cm]

global White_space Gray_space Black_space
White_space = [];
Gray_space = [];
Black_space = [];
White_space = ones(1,L_street_binning); % percent
Gray_space = zeros(1,L_street_binning); % percent
Black_space = zeros(1,L_street_binning); % percent

park_que = 0;
parking_cars = 0;
time = 0;

full_or_part = false; % 0-full 1-partial
with_marks = true;
full_marks = false; % 0-half street 1-full street
to_plot_illustration = false;
L_car_to_promote = 4.5; % in meters

if full_or_part %P
    if with_marks %M
        if full_marks %F
            case_study = 'PMF';
        else
            case_study = 'PM';
        end
    else %nM
        case_study = 'PnM';
    end
else %F
    if with_marks %M
        if full_marks %F
            case_study = 'FMF';
        else
            case_study = 'FM';
        end
    else %nM
        case_study = 'FnM';
    end
end

% Define car
L_car_width = 2;
L_car_stearing = 1.15; % [percent] additional needed space

if full_or_part
    L_car_verity = [7 4 3 2.5 2 1.5]; % [meters]
    P_car_disterbution = [0.02 0.2 0.18 0.3 0.28 0.02]; % disterbution out of ten cars
    car_min_space = min(L_car_verity .* L_car_stearing);
else
    L_car_verity = [1:0.2:11]; % [meters]
    P_car_disterbution1 = [length(L_car_verity):-1:1]./sum(length(L_car_verity):-1:1); % disterbution out of ten cars
    P_car_disterbution2 = [1:length(L_car_verity)]./sum(1:length(L_car_verity)); % disterbution out of ten cars
    P_car_disterbution3 = ones(1,length(L_car_verity))./length(L_car_verity); % disterbution out of ten cars
    P_car_disterbution = P_car_disterbution3;
    car_min_space = sum(P_car_disterbution.* L_car_stearing .* L_car_verity);
end

marked_park_size = round(car_min_space*2.0 * 100) * with_marks;
marked_park_size = round(car_min_space*20.0 * 100) * with_marks;
marked_park_size = round(L_car_to_promote*100) * with_marks;

disp(['The marks size in meters is : ' num2str(marked_park_size/100)])
if with_marks
    marked_park_size_index = 1;
    if full_marks
        disp(['The street is fully marked'])
        while marked_park_size_index<length(Black_space)
            Black_space(marked_park_size_index) = 1;
            marked_park_size_index = marked_park_size_index + marked_park_size;
        end

    else
        disp(['The street is half marked'])
        while marked_park_size_index<length(Black_space)/2
            Black_space(marked_park_size_index) = 1;
            marked_park_size_index = marked_park_size_index + marked_park_size;
        end
    end
    White_space = White_space - Black_space;
end
% driver defined in class
N_drivers = 5000; % optimal 5000
L_car = randsample(L_car_verity,N_drivers,true,P_car_disterbution);
L_car_min_space = L_car .* L_car_stearing;

for i=1:N_drivers
    driver(i) = Driver_obj(i);
end

%% calculation
% plotting street
if to_plot_illustration
    figure(1)
    clf
    if sum(Black_space)
        x_vec = find(Black_space~=0);
        for i=1:length(x_vec)
            line([x_vec(i) x_vec(i)]./100, [2 -0.3])
        end
        axis([-0.1 L_street+0.1 -5 10])
        title({'Illustration of the Steet',['time = ' num2str(time)]})
        xlabel('Place along the street [m]')
    else
        axis([-0.1 L_street+0.1 -5 10])
        title({'Illustration of the Steet',['time = ' num2str(time)]})
        xlabel('Place along the street [m]')
    end
    patch([0 L_street L_street 0],[L_car_width*3 L_car_width*3 L_car_width L_car_width]+0.2,[.1 .1 .1]);
    street_marks=0;
    while street_marks<L_street
        patch([street_marks street_marks+0.5 street_marks+0.5 street_marks], ...
            [L_car_width*2+0.1 L_car_width*2+0.1 L_car_width*2-0.1 L_car_width*2-0.1]+0.2,[1 1 1]);
        street_marks = street_marks + 1.5;
    end
end

cars_still_park = true;
d_id = 1;
while cars_still_park
    % new driver
    driver(d_id).parked_position = finelyzing_parking(find_park(L_car_min_space(d_id) + L_car_min_space(d_id) * driver(d_id).additional_space), L_car(d_id), L_car_min_space(d_id) + L_car_min_space(d_id) * driver(d_id).additional_space);
    
    time = time + driver(d_id).arrival_time;
    driver(d_id).abs_arrival_time = time;
    
    if driver(d_id).parked_position < 0
        gemming_driver(d_id) = d_id;
        gemming_cars(d_id) = L_car(d_id);
        gemming_time(d_id) = driver(d_id).waiting_time;
        
        driver(d_id).waiting_que = true;
        
        driver(d_id).abs_leaving_time = driver(d_id).waiting_time + time;
    else
        driver(d_id).parked_position_end = driver(d_id).parked_position + L_car(d_id)*100;
        
        White_space(int16(driver(d_id).parked_position:driver(d_id).parked_position_end)) = 0;
        Gray_space(int16(driver(d_id).parked_position:driver(d_id).parked_position_end)) = 1;
        driver(d_id).parking_bool = true;
        
        driver(d_id).abs_leaving_time = driver(d_id).leaving_time + time;
    end
    
    % plotting updated street
    if to_plot_illustration && driver(d_id).parking_bool
        car_p(d_id) = patch([driver(d_id).parked_position driver(d_id).parked_position_end ...
            driver(d_id).parked_position_end driver(d_id).parked_position]./100,[L_car_width L_car_width 0 0],[.5 .5 .5]);
        if time<60
            title({'Illustration of the Steet',['time = ' num2str(time) ' min']})
        elseif time<60*24
            title({'Illustration of the Steet',['time = ' num2str(time/60) ' hours']})
        elseif time<60*24*7
            title({'Illustration of the Steet',['time = ' num2str(time/60/24) ' days']})
        elseif time<60*24*7*4
            title({'Illustration of the Steet',['time = ' num2str(time/60/24/7) ' weeks']})
        elseif time<60*24*365
            title({'Illustration of the Steet',['time = ' num2str(time/60/24/30) ' months']})
        end
    end
    
    %leaving driver
    if ~isempty(find([driver(1:d_id).abs_leaving_time] < time))
        index_left_street = [driver(([driver(1:d_id).abs_leaving_time] < time).*[driver(1:d_id).parking_bool] == true).driver_id];
        index_left_street1 = [driver(([driver(1:d_id).abs_leaving_time] < time).*[driver(1:d_id).waiting_que] == true).driver_id];
        if to_plot_illustration
            delete(car_p(index_left_street))
        end
        for i=[index_left_street index_left_street1]
            if driver(i).parking_bool
                White_space(int16(driver(i).parked_position:driver(i).parked_position_end)) = 1;
                Gray_space(int16(driver(i).parked_position:driver(i).parked_position_end)) = 0;
                driver(i).parking_bool = false;
            end
            if driver(i).waiting_que
                driver(i).waiting_que = false;
            end
        end
    end
    
    parking_cars(d_id) = sum([driver(:).parking_bool]);
    park_que(d_id) = sum([driver.waiting_que]);
    time_vec(d_id) = time;
    if to_plot_illustration
        pause(0.05)
    else
        if mod(d_id,100) == 0
            waitbar(d_id/N_drivers)
        end
    end
    d_id = d_id+1;
    if d_id>=N_drivers+1
        break
    end
end

[uniques_L,numUnique_L] = count_unique(L_car);
gemming_cars(gemming_cars==0)=[];
[uniques_G,numUnique_G] = count_unique(gemming_cars);

numUnique_G = [numUnique_G ; zeros(length(setdiff(uniques_L,uniques_G)),1)];
uniques_G = [uniques_G ; setdiff(uniques_L,uniques_G)];

[uniques_G, uniques_G_order] = sort(uniques_G);
numUnique_G = numUnique_G(uniques_G_order,:);

missed_park_by_lenght_percent = numUnique_G./numUnique_L.*100;

%% Numerical
% sigmoid function (“S” shape curve)
a = 100;

% for full lenght spectrum no marks
b = 1.474; % (1.389, 1.559)
c = 17.51; % (14.35, 20.67)
P_to_park_com = @(car_length) a ./ (1 + exp(-1 .* car_length ./b )) .^c;

% for partial lenght spectrum with marks
b = 0.4825; % (0.06834, 0.8966)
c = 192.6; % (-744.4, 1130)
P_to_park = @(car_length) a ./ (1 + exp(-1 .* car_length ./b )) .^c;

% Define Start points, fit-function and fit curve
x0 = [b c];
fitfun = fittype( @(b,c,x) a./(1 + exp(-1.*x./b)).^c );
[fitted_curve,gof] = fit(uniques_L,missed_park_by_lenght_percent,fitfun,'StartPoint',x0);
% Save the coeffiecient values for a,b,c and d in a vector
coeffvals = coeffvalues(fitted_curve);


%% results
gemming_time(gemming_cars==0)=[];
gemming_cars(gemming_cars==0)=[];
gemming_cars_num = cumsum(gemming_cars~=0);

fit_ps_ot = fit(time_vec(100:end)',parking_cars(100:end)','poly1');
fit_ps_ot.p1 = 0;

fit_wc_ot = fit(time_vec(100:end)',park_que(100:end)','poly1');
fit_wc_ot.p1 = 0;

SEM = std(parking_cars)/sqrt(N_drivers);

figure(2)
bar(uniques_L,missed_park_by_lenght_percent)
hold on
line([marked_park_size marked_park_size]./100,[0 100],'Color','red')
grid on
title('Probability of Cars to Missed Parking by Length')
xlabel('Length [m]')
ylabel('Percent of Cars Missed Parking')
hold off
legend({'Cars missed Parking','Parking spots size'},'Location','east')
eval(['print -depsc P2P' case_study])

figure(3)
plot(1:0.1:max(uniques_L) , fitted_curve(1:0.1:max(uniques_L)),'m-x')
hold on
line([marked_park_size marked_park_size]./100,[0 100],'Color','red')
plot(1:0.1:max(uniques_L) , P_to_park_com(1:0.1:max(uniques_L)),'b-+')
% plot(uniques_L , P_to_park(uniques_L),'ko','markerfacecolor','k')
% plot(1:0.1:max(uniques_L) , P_to_park(1:0.1:max(uniques_L)),'g-s')
title('Probability of Cars to Missed Parking by Length')
xlabel('Length [m]')
ylabel('Percent of Cars Missed Parking')
hold off
legend({'Parking probability for this case','Parking spots size','Parking probability for full spectrum case','Fitted probability for parking this case','Fitted probability for parking without marks'},'Location','southeast')
eval(['print -depsc P2Pfun' case_study])

figure(4)
[y,x] = hist(L_car,uniques_L); % total
[y1,x1] = hist(L_car([driver(:).parked_position]>0),uniques_L); % total parked
[y2,x3] = hist(gemming_cars,uniques_L); % total missed
bar3(x,[y2;y1;y]')
title('Number of Cars Vs Length')
xticklabels({'Missed Parking Popolation','Parking Popolation','Total Popolation'})
ylabel('Length [m]')
set(get(gca(),'Children'),'FaceAlpha',0.4)
view(-100,30)
eval(['print -depsc Dis' case_study])

figure(5)
plot(time_vec,parking_cars,'.')
hold on
plot(time_vec,park_que,'r.')
plot(fit_ps_ot,'b')
plot(fit_wc_ot,'r')
hold off
legend({'Parking Cars','Flux of Cars With no Parking [car/min]',['Mean Number Of Parking Cars is ' num2str(fit_ps_ot.p2)],['Mean Number Of Waiting Cars is ' num2str(fit_wc_ot.p2)]},'Location','west')
grid on
title({'Number of Parking Cars vs. Time'})
xlabel('Time [min]')
ylabel('Number of Cars')
eval(['print -depsc r' case_study])

if 0
    %figure(6)
    openfig('street.fig');
    print -depsc street
    
    figure(7)
    hist([driver.arrival_time],100)
    title('Arrival Time Distribution')
    xlabel('Time [min]')
    ylabel('Number of Cars')
    print -deps arrival
    
    figure(8)
    hist([driver.leaving_time],100)
    title('Leaving Time Distribution')
    xlabel('Time [min]')
    ylabel('Number of Cars')
    print -deps leaving
    
    figure(9)
    hist([driver.waiting_time],100)
    title('Waiting Time Distribution (in case of no parking)')
    xlabel('Time [min]')
    ylabel('Number of Cars')
    print -deps waiting
    
    figure(10)
    hist([driver.additional_space],100)
    title('Additional Space Needed By Driver')
    xlabel('Additional Needed Space [meter]')
    ylabel('Number of Cars')
    print -deps additional
    
    if full_or_part
        figure(11)
        bar(L_car_verity,P_car_disterbution.*100)
        title('Car Disterbution By Length')
        axis([0 max(L_car_verity)+1 0 1.5*max(P_car_disterbution.*100)])
        xlabel('Length [meter]')
        ylabel('Percentage [%]')
        print -deps car_l
    else
        figure(11)
        bar(L_car_verity,P_car_disterbution1.*100)
        title('Car Disterbution By Length')
        axis([0 max(L_car_verity)+1 0 1.5*max(P_car_disterbution1.*100)])
        xlabel('Length [meter]')
        ylabel('Percentage [%]')
        print -deps car_l_f1
        
        figure(12)
        bar(L_car_verity,P_car_disterbution2.*100)
        title('Car Disterbution By Length')
        axis([0 max(L_car_verity)+1 0 1.5*max(P_car_disterbution2.*100)])
        xlabel('Length [meter]')
        ylabel('Percentage [%]')
        print -deps car_l_f2
        
        figure(13)
        bar(L_car_verity,P_car_disterbution3.*100)
        title('Car Disterbution By Length')
        axis([0 max(L_car_verity)+1 0 1.5*max(P_car_disterbution3.*100)])
        xlabel('Length [meter]')
        ylabel('Percentage [%]')
        print -deps car_l_f3
    end
end

disp(['The mean value of parking cars is ' num2str(mean(mean(parking_cars))) ' cars.'])
disp(['The Standard Error Of The Mean is ' num2str(SEM)])
disp(' ')
disp(['E(parking cars)=(' num2str(mean(parking_cars)-SEM) ',' num2str(mean(parking_cars)+SEM) ')'])

disp(['The street max capacity was ' num2str(max(parking_cars)) ' parking cars.'])
disp(['The street min capacity was ' num2str(min(parking_cars(100:end))) ' parking cars.'])
disp(['Street of ' num2str(L_street) ' meters ' num2str(mean(parking_cars)) ' parking cars.'])
disp(['Free capacity is ' num2str(sum(White_space)/(L_street*100)*100) '% on average ' num2str(floor(sum(White_space)/(car_min_space*100))) ' more cars could park.'])
disp(['Ocupide capacity on average ' num2str(mean(parking_cars)/max(parking_cars)*100) '.'])
disp(['Blocked capacity is ' num2str(sum(Black_space)/(L_street*100)*100) '%.'])

disp 'Done ...'

%% functions
function car_parked = find_park(needed_space_to_park)
global White_space Gray_space Black_space
needed_space_in_bins = round(needed_space_to_park * 100);

ES = movsum(White_space,[0 needed_space_in_bins-1],'Endpoints','discard');
PS = find(ES==needed_space_in_bins);
if isempty(PS)
    car_parked=-1e5;
    return
else
    car_parked = PS(randi(numel(PS)));
    car_parked =  car_parked/100;
    return
end
end

function car_parked = finelyzing_parking(space_to_park, car_size, avail_space)
if space_to_park<0
    car_parked=-1e5;
    return
end

car_size = car_size * 100;
avail_space = floor(avail_space*100);

car_parked = space_to_park * 100 + randi(int16(avail_space-car_size));
end