function new_param_axis = create_axis_ticks(param, param_ticks_steps)

	param = round(param, 2);
	param_str = arrayfun(@num2str, param, 'UniformOutput', false);

	param_axis_indices = 1:param_ticks_steps:length(param);
	param_axis = param_str(param_axis_indices);

	if param_ticks_steps-1 == 0
		new_param_axis = cell(1, length(param));
		for i = 1:length(param)
			idx = (i-1)*param_ticks_steps + 1;
			new_param_axis(i) = param_axis(i);
		end
	else
		n_y_axis_empty_cells = param_ticks_steps-1;
		new_param_axis = cell(1, length(param));
		y_axis_empty_cells = repmat({' '}, 1, n_y_axis_empty_cells);
		for i = 1:param_ticks_steps
			idx = (i-1)*param_ticks_steps + 1;
			new_param_axis(idx:idx+param_ticks_steps-1) = ...
				[param_axis(i), y_axis_empty_cells];
		end
	end

end
