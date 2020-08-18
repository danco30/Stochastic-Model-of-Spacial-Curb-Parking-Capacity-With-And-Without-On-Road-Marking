% define driver
classdef Driver_obj
    properties
        driver_id
        arrival_time
        abs_arrival_time
        leaving_time
        abs_leaving_time
        waiting_time
        additional_space
        preffered_spot % not used
        parked_position
        parked_position_end
        to_park_time % for traffic delays not used
        waiting_que = false;
        parking_bool = false;
    end
    methods
        function obj = Driver_obj(id)
            obj.driver_id = id;
            obj.arrival_time = exprnd(15,1); % [min] since last
            obj.leaving_time = obj.erlangrnd(); % since arraival
            obj.waiting_time = exprnd(50,1); % [min] assume FIFO
            obj.additional_space = 0.05*rand(1); % [percent of car size] due to conidence
            obj.to_park_time = exprnd(5,1); % [min]
            obj.preffered_spot =  mean(rand(1,2)); % in the street in percentadge
        end
    end
    methods(Static)
        function rnd_x = erlangrnd()
            mu = [120 200];
            k = [2 10];
            
            short_p = 0.6;
            long_p = 0.4;
            
            short_erlang_fun1 = @(x) x.^(k(1)-1).*exp(-x/mu(1))/(mu(1).^k(1).*factorial(k(1)-1)).*short_p;
            long_erlang_fun2 = @(x) x.^(k(2)-1).*exp(-x/mu(2))/(mu(2).^k(2).*factorial(k(2)-1)).*long_p;
            general_erlang_fun =@(x) short_erlang_fun1(x) + long_erlang_fun2(x);
            
            rnd_x = randpdf(general_erlang_fun(0:5000),0:5000,[1 1]);
        end
    end
end