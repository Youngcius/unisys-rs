pub fn is_power_of_two(num: usize) -> bool {
    (num & (num - 1) == 0) && num != 0
}

