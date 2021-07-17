import random


def get_instructions(reference_size, repetitive_region_size, coverage, seed = 49):
    num_rr = int(reference_size * coverage / repetitive_region_size)
    DISTANCE = 10
    random.seed(seed)
    start_positions = []
    while len(start_positions) < num_rr:
        position = random.randint(0, reference_size - 1 - repetitive_region_size)
        include = True
        for start_position in start_positions:
            if abs(position - start_position) <= repetitive_region_size + DISTANCE:
                include = False
                break
        if include:
            start_positions.append(position)

    instructions = {}
    for i, start_position in enumerate(start_positions):
        end_position = start_position + repetitive_region_size
        instructions[start_position] = seed
        instructions[end_position] = seed + (i + 1) * seed
    return instructions
