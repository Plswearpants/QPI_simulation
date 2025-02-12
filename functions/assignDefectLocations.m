function defect_locations = assignDefectLocations(num_defects, lattice_size)
    % Randomly assign defect locations within the lattice without overlapping
    defect_locations = zeros(num_defects, 2);
    assigned_positions = false(lattice_size, lattice_size); % Logical array to track assigned positions
    
    for i = 1:num_defects
        while true
            new_location = randi([1, lattice_size], 1, 2);
            if ~assigned_positions(new_location(1), new_location(2))
                defect_locations(i, :) = new_location;
                assigned_positions(new_location(1), new_location(2)) = true;
                break;
            end
        end
    end
end
