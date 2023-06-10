function param_str = param2str(param);
	param_str = num2str(param);
	% only if coupling is not zero, we use letter sequence as of third letter in coupling_string
	if param ~= 0
		param_str = param_str(3:end);
	end
			
	% if coupling_string is only one number as string, add zero
	if length(param_str) == 1;
		param_str = [param_str, '0'];
	end
end