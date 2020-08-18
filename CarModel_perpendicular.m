% curb parking
% three kinds of spaces:
% white - available
% gray - parking
% black - obstile
clear
clc
clf

format longg
rng(123); % creating a single seed for repetetivity

% Define street
L_street = 10; % [meters]
L_street_binning = L_street*100; % [cm]

global White_space Gray_space Black_space
White_space = [];
Gray_space = [];
Black_space = [];
White_space = ones(1,L_street_binning); % percent
Gray_space = zeros(1,L_street_binning); % percent
Black_space = zeros(1,L_street_binning); % percent

waiting_cars = 0;
parking_cars = 0;
missed_park = 0;

time = 0;
with_marks = false;
to_plot_illustration = true;

% Define car
L_car_width = 4;
L_car_length_verity = [2 2.5]; % length [meters]
L_car_stearing = 1.15; % [percent] additional needed space
P_car_disterbution = [0.8 0.2]; % disterbution out of ten cars

car_min_space = sum(P_car_disterbution.* L_car_stearing .* L_car_length_verity);
if with_marks
    marked_park_size = round(car_min_space*1.23 * 100);
    marked_park_size_index = 1;
    while marked_park_size_index<length(Black_space)
        Black_space(marked_park_size_index) = 1;
        marked_park_size_index = marked_park_size_index + marked_park_size + randi(50);
    end
    White_space = White_space - Black_space;
end
% driver defined in class
N_drivers = 300;
L_car = randsample(L_car_length_verity,N_drivers,true,P_car_disterbution);
L_car_min_space = L_car .* L_car_stearing;

for i=1:N_drivers
    driver(i) = Driver_obj(i);
end

%% calculation
% plotting street
if to_plot_illustration
    figure(1)
    if sum(Black_space)
        x_vec = find(Black_space~=0);
        for i=1:length(x_vec)
            line([x_vec(i) x_vec(i)]./100, [2 -2])
        end
        axis([-0.1 L_street+0.1 -5 10])
        title({'Illustration of the Steet',['time = ' num2str(time)]})
        xlabel('Place along the street [m]')
    else
        axis([-0.1 L_street+0.1 -5 10])
        title({'Illustration of the Steet',['time = ' num2str(time)]})
        xlabel('Place along the street [m]')
    end
end
cars_still_park = true;
d_id = 1;
while cars_still_park
    if d_id>=N_drivers
        break
    end
    
    % new driver
    driver(d_id).parked_position = finelyzing_parking(find_park(L_car_min_space(d_id) + L_car_min_space(d_id) * driver(d_id).additional_space), L_car(d_id), L_car_min_space(d_id) + L_car_min_space(d_id) * driver(d_id).additional_space);
    driver(d_id).parked_position_end = driver(d_id).parked_position + L_car(d_id)*100;
    
    time = time + driver(d_id).arrival_time;
    driver(d_id).abs_arrival_time = time;
    driver(d_id).abs_leaving_time = driver(d_id).leaving_time + time;
    
    if driver(d_id).parked_position<0
        time = time + driver(d_id).waiting_time;
        gemming_cars(d_id) = L_car(d_id);
        gemming_time(d_id) = time;
    else
        White_space(int16(driver(d_id).parked_position):int16(driver(d_id).parked_position_end)) = 0;
        Gray_space(int16(driver(d_id).parked_position):int16(driver(d_id).parked_position_end)) = 1;
        driver(d_id).parking_bool = true;
    end
    
    % plotting updated street
    if to_plot_illustration
        car_p(d_id) = patch([driver(d_id).parked_position driver(d_id).parked_position_end ...
            driver(d_id).parked_position_end driver(d_id).parked_position]./100,[L_car_width L_car_width 0 0],[.5 .5 .5]);
        title({'Illustration of the Steet',['time = ' num2str(time)]})
    end
    
    %leaving driver
    if ~isempty(find([driver(1:d_id).abs_leaving_time] < time))
        index_left_street = [driver([driver(1:d_id).abs_leaving_time] < time).driver_id];
        if to_plot_illustration
            delete(car_p(index_left_street))
        end
        for i=index_left_street
            if driver(i).parking_bool
                White_space(int16(driver(i).parked_position):int16(driver(i).parked_position_end)) = 1;
                Gray_space(int16(driver(i).parked_position):int16(driver(i).parked_position_end)) = 0;
                driver(i).parking_bool = false;
            end
        end
    end
    
    d_id = d_id+1;
    parking_cars(d_id) = sum([driver(:).parking_bool]);
    time_vec(d_id) = time;
    pause(0.05)
end

%% results
gemming_time(gemming_cars==0)=[];
gemming_cars(gemming_cars==0)=[];
gemming_cars_num = cumsum(gemming_cars~=0);

fit_ps_ot = fit(time_vec(100:end)',parking_cars(100:end)','poly1');

SEM = std(parking_cars)/sqrt(N_drivers);
figure(2)
hist(gemming_cars)
grid on
title('Number of Parking Cars vs. Time')
xlabel('Time [min]')
ylabel('Number of Cars')

figure(3)
plot(time_vec,parking_cars,'.')
hold on
plot(gemming_time,[0 diff(gemming_cars_num)],'r.')
plot(fit_ps_ot)
hold off
legend('Parking Cars','Flux of Cars With no Parking [car/min]','Fitting the stedy state')

grid on
title({'Number of Parking Cars vs. Time',['Degration per min = ' num2str(fit_ps_ot.p1)]})
xlabel('Time [min]')
ylabel('Number of Cars')

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
    car_parked=-1;
    return
else
    car_parked = PS(randi(numel(PS)));
    car_parked =  car_parked/100;
    return
end
end

function car_parked = finelyzing_parking(space_to_park, car_size, avail_space)
if space_to_park==0
    car_parked=0;
    return
end

car_size = car_size * 100;
avail_space = floor(avail_space*100);

car_parked = space_to_park * 100 + randi(avail_space-car_size);
end