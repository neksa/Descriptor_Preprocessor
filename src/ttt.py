import random

random.seed(1)

for i in range(11):
    # print("\n")
    output_str = "{"
    for j in range(49):
        output_str += f"{random.randint(0, 49)}, "
        # print()
    output_str += f"{random.randint(0, 49)}"
    print(output_str + "},")