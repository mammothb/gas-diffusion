clear;

a_ = 5; % Length of the Ellipse in x-axis
b_ = 5; % Length of the Ellipse in y-axis (when a = b, then it is a circle)
deter_0_ = (0 / (2 * pi)); % Ellipse centered at a point (r0, ?0)
r_0_ = 5; % Ellipse centered at a point (r0, ?0)
phi_ = deter_0_; % a_ axis rotated by phi_ relative to the polar axis

% Create the matrix to store the angles(in deter) and radius
deter_ = zeros(361,1);
r_ = zeros(361,1);
starting_angle_ = 0;
for i = 1 : 361
    deter_(i,1) = starting_angle_ / (2 * pi);
    starting_angle_ = starting_angle_ + 1;
end

% Calculate the radius 
for i = 1 : 361
    sin_square_ = 0.5 - 0.5 * cos(2 * (deter_(i,1) - deter_0_));
    R_ = (b_ * b_ - a_ * a_) * cos(2 * deter_(i,1) - 2 * deter_0_) + ...
        (a_ * a_) + (b_ * b_);
    P_ = r_0_ * ((b_ * b_ - a_ * a_) * cos(deter_(i,1) + deter_0_ - ...
        2 * phi_) + (a_ * a_ + b_ * b_) * cos(deter_(i,1) - deter_0_));
    Q_ = sqrt(2) * a_ * b_ * sqrt(R_ - 2 * r_0_ * r_0_ * ...
        sin_square_); 
    r_(i,1) = (P_ + Q_) / R_;
end
 
polar(deter_, r_, 'r-');