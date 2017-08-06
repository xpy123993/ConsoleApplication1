#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>

#include <algorithm>

#define VALUE_DEFAULT -1

using namespace std;

//Parameters for test data
const int ALLOCATE_MAX_M = 10;
const int ALLOCATE_MAX_N = 400;

const int JOB_PROCTIME_MIN = 1;
const int JOB_PROCTIME_MAX = 99;

const char* JOB_FILE = "test.dat";

//Parameters for batch test
const bool GENERATE_JOB_ON_START = true;
const bool PRINT_ITERATION_SOLUTION = false;
const bool PRINT_TABLE = true;

const int Ms[] = { 3, 5, 10 };
const int NdMs[] = { 5, 10, 20 };
const int Bs[] = { 20, 30, 50 };
const int ITERATION_IN_EACH_CASE = 1;

//Parameters for pso
const double MUTATION_PROBABILITY = 0.02;
const double SELECTION_PROBABILITY = 0.5;
const double SAVE_PROBABILITY = 0.40;

const int POPULATION_SIZE = 50;
const int PSO_ITERATION_MAXN = 10000;
const int PSO_UNCHANGED_ITERATION_MAXN = 2000;
const bool PSO_UNCHANGED_ITERATION_MODE = true;

//Parameter for tabu search
const int TABU_LIST_LENGTH = 20;
const int TABU_ITERATION_MAXN = 10000;
const int TABU_UNCHANGED_ITERATION_MAXN = 2000;
const bool TABU_UNCHANGED_ITERATION_MODE = true;
const int TABU_RANDOM_SELECTION_ITERATION = 50;

const int NUM_MAXN = 1000007;

int M = 5, N = 100, b = 20;
double r = 0.5, t = 0.6;

struct Job
{
	int proc_time[ALLOCATE_MAX_M];
	int due_date;
};

Job test_data[ALLOCATE_MAX_N];

int indexof(int* array, int size, int value)
{
	for (int i = 0; i < size; i++)
	{
		if (array[i] == value)
			return i;
	}
	return -1;
}

int argmax(int* array, int size)
{
	int max_index = 0;
	for (int i = 1; i < size; i++)
	{
		if (array[i] > array[max_index])
			max_index = i;
	}
	return max_index;
}

int argmin(int* array, int size)
{
	int min_index = 0;
	for (int i = 1; i < size; i++)
	{
		if (array[min_index] > array[i])
			min_index = i;
	}
	return min_index;
}

double avg(int* array, int size)
{
	double sum = 0;
	for (int i = 0; i < size; i++)
		sum += array[i];
	return sum / size;
}

int cmp_edd(int index_l, int index_r)
{
	return test_data[index_l].due_date < test_data[index_r].due_date;
}

int cmp_avg_proc(const int index_l, const int index_r)
{
	int sum = 0;
	for (int i = 0; i < M; i++)
		sum += test_data[index_l].proc_time[i] - test_data[index_r].proc_time[i];
	return sum < 0;
}

void draw_line(int line_width)
{
	while (line_width--)
		printf("-");
	printf("\n");
}

void print_job_info(Job& j)
{
	draw_line(10);
	printf("Array[] Process Time:\n");
	for (int i = 0; i < M; i++)
	{
		printf("%d ", j.proc_time[i]);
	}
	printf("\n");
	printf("Due Date: %d\n", j.due_date);
	draw_line(10);
	printf("\n");
}

int u(int min, int max)
{
	return min + (rand() % (max - min));
}

void store_test_data()
{
	FILE *fp = fopen(JOB_FILE, "wb");
	if (!fp) return;

	fwrite(test_data, sizeof(Job), N, fp);
	fclose(fp);
}

void load_test_data()
{
	FILE *fp = fopen(JOB_FILE, "rb");
	if (!fp) return;

	fread(test_data, sizeof(Job), N, fp);
	fclose(fp);
}

/*
2017-5-9
数据生成方法
参考论文
Two machine flow shop scheduling to minimize total late work
*/
void generate_test_data_2()
{
	int mapper_index_1, mapper_index_2;
	int mapper[ALLOCATE_MAX_N];
	double min_due_date, max_due_date;
	for (int i = 0; i < N; i++)
	{
		mapper[i] = i;
		for (int j = 0; j < M; j++)
		{
			test_data[i].proc_time[j] = u(JOB_PROCTIME_MIN, JOB_PROCTIME_MAX);
		}
	}

	sort(mapper, mapper + N, cmp_avg_proc);

	for (int i = 0; i < N; i++)
	{
		mapper_index_1 = mapper[i];
		min_due_date = avg(test_data[mapper_index_1].proc_time, M);
		max_due_date = min_due_date;

		for (int k = 0; k < i; k++)
		{
			mapper_index_2 = mapper[k];
			max_due_date += avg(test_data[mapper_index_2].proc_time, M) / b;
		}

		test_data[mapper_index_1].due_date = u(min_due_date, max_due_date);
	}
	store_test_data();
}

/*
2017-2-27
生成数据参考论文
Parallel-machine scheduling to minimize tardiness penalty and power cost
*/
void generate_test_data()
{

	generate_test_data_2();
	return;


	int proc_time_sum = 0, min_due_date, max_due_date, due_date_range;

	double avg_proc_time;
	double avg_machine_load = 0;

	for (int i = 0; i < N; i++)
	{
		//for every test-job
		proc_time_sum = 0;
		for (int j = 0; j < M; j++)
		{
			// generate individual speed for every machine
			test_data[i].proc_time[j] = u(JOB_PROCTIME_MIN, JOB_PROCTIME_MAX);
			proc_time_sum += test_data[i].proc_time[j];
		}
		avg_proc_time = 1.0 * proc_time_sum / M;
		avg_machine_load += avg_proc_time;
	}

	avg_machine_load = avg_machine_load / M;

	due_date_range = avg_machine_load * r;

	min_due_date = avg_machine_load * (1.0 - t - r / 2);
	max_due_date = min_due_date + due_date_range;

	for (int i = 0; i < N; i++)
	{
		// generate due date
		test_data[i].due_date = u(min_due_date, max_due_date);
		//print_job_info(test_data[i]);
	}

	//printf("due date: [%d, %d)\n", min_due_date, max_due_date);

	store_test_data();
}

struct Solution
{
public:
	int value[ALLOCATE_MAX_M][ALLOCATE_MAX_N];
	Solution()
	{
		for (int i = 0; i < ALLOCATE_MAX_M; i++)
			for (int j = 0; j < ALLOCATE_MAX_N; j++)
				value[i][j] = VALUE_DEFAULT;
	}

	bool equals(const Solution& solution)
	{
		for (int i = 0; i < M; i++)
		{
			if (memcmp(value[i], solution.value[i], N * sizeof(int)) != 0)
				return false;
		}
		return true;
	}
};

int get_lateness(int job_index, int machine_index, int completion_time)
{
	return min(test_data[job_index].proc_time[machine_index], max(0, completion_time - test_data[job_index].due_date));
}

int evaluate(const Solution& solution)
{
	int proc_time[ALLOCATE_MAX_M] = { 0 };
	int latework = 0;

	int job_detected_count = 0;

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (~solution.value[i][j])
			{
				Job& current_job = test_data[solution.value[i][j]];
				proc_time[i] += current_job.proc_time[i];
				latework += get_lateness(solution.value[i][j], i, proc_time[i]);
				job_detected_count++;
			}
		}
	}

	if (job_detected_count != N)
	{
		printf("ERROR: Job data lost, current number: %d (%d)\n", job_detected_count, N);
	}

	return latework;
}

bool check_solution(const Solution& solution)
{
	int jobs[ALLOCATE_MAX_N];
	int count = 0;
	for (int i = 0; i < N; i++)
	{
		jobs[i] = -1;
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (solution.value[i][j] != VALUE_DEFAULT)
			{
				count++;
				jobs[solution.value[i][j]] ++;
				if (jobs[solution.value[i][j]] > 1) return false;
			}
		}
	}

	return count == N;
}

void print_solution(const Solution& solution, const char* tagName = NULL)
{
	draw_line(10);
	if (tagName != NULL)
	{
		printf("Solution: %s\n", tagName);
	}
	int proc_time[ALLOCATE_MAX_M] = { 0 };
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (solution.value[i][j] != VALUE_DEFAULT)
				proc_time[i] += test_data[solution.value[i][j]].proc_time[i];
		}
	}
	for (int i = 0; i < M; i++)
	{
		printf("M%d (Process Time = %d): ", i, proc_time[i]);
		for (int j = 0; j < N; j++)
		{
			if (solution.value[i][j] != VALUE_DEFAULT)
			{
				Job &current_job = test_data[solution.value[i][j]];
				printf("%d ", solution.value[i][j]);
			}
		}
		printf("\n");

	}

	printf("Value = %d\n", evaluate(solution));
	draw_line(10);
}



class ListOptimizer
{
public:
	Solution get_min_makespan_first_solution()
	{
		Solution solution;
		int proc_time[ALLOCATE_MAX_M] = { 0 };
		int proc_count[ALLOCATE_MAX_M] = { 0 };
		int selected_index = 0;
		for (int i = 0; i < N; i++)
		{
			selected_index = argmin(proc_time, M);
			proc_time[selected_index] += test_data[i].proc_time[selected_index];
			solution.value[selected_index][proc_count[selected_index]++] = i;
		}
		return solution;
	}
	Solution get_min_proc_time_first_solution()
	{
		Solution solution;
		int proc_count[ALLOCATE_MAX_M] = { 0 };
		int selected_index = 0;
		for (int i = 0; i < N; i++)
		{
			selected_index = argmin(test_data[i].proc_time, M);
			solution.value[selected_index][proc_count[selected_index]++] = i;
		}
		return solution;
	}
	Solution get_minc_solution()
	{
		Solution solution;
		int proc_time[ALLOCATE_MAX_M] = { 0 };
		int proc_count[ALLOCATE_MAX_M] = { 0 };
		int selected_index = 0;
		for (int i = 0; i < N; i++)
		{
			selected_index = argmin(proc_time, M);
			proc_time[selected_index] += test_data[i].proc_time[selected_index];
			solution.value[selected_index][proc_count[selected_index]++] = i;
		}
		if (!check_solution(solution))
		{
			printf("Warning: Unusual Solution found\n");
			print_solution(solution);
		}
		return solution;
	}
	Solution get_miny_solution()
	{
		Solution solution;
		int proc_time[ALLOCATE_MAX_M] = { 0 };
		int proc_count[ALLOCATE_MAX_M] = { 0 };
		int proc_late[ALLOCATE_MAX_M] = { 0 };
		int selected_index = 0;
		for (int i = 0; i < N; i++)
		{
			selected_index = argmin(proc_late, M);
			proc_time[selected_index] += test_data[i].proc_time[selected_index];
			proc_late[selected_index] += get_lateness(i, selected_index, proc_time[selected_index]);
			solution.value[selected_index][proc_count[selected_index]++] = i;
		}
		if (!check_solution(solution))
		{
			printf("Warning: Unusual Solution found\n");
			print_solution(solution);
		}
		return solution;
	}
	Solution get_spt_edd_solution()
	{
		Solution solution;
		int proc_count[ALLOCATE_MAX_M] = { 0 };
		int selected_index;
		for (int i = 0; i < N; i++)
		{
			selected_index = argmin(test_data[i].proc_time, M);
			solution.value[selected_index][proc_count[selected_index]++] = i;
		}
		for (int i = 0; i < M; i++)
		{
			if (proc_count[i] > 0)
				sort(solution.value[i], solution.value[i] + proc_count[i], cmp_edd);
		}
		if (!check_solution(solution))
		{
			printf("Warning: Unusual Solution found\n");
			print_solution(solution);
		}
		return solution;
	}
	Solution get_lpt_edd_solution()
	{
		Solution solution;
		int proc_count[ALLOCATE_MAX_M] = { 0 };
		int selected_index;
		for (int i = 0; i < N; i++)
		{
			selected_index = argmax(test_data[i].proc_time, M);
			solution.value[selected_index][proc_count[selected_index]++] = i;
		}
		for (int i = 0; i < M; i++)
		{
			if (proc_count[i] > 0)
				sort(solution.value[i], solution.value[i] + proc_count[i], cmp_edd);
		}
		if (!check_solution(solution))
		{
			printf("Warning: Unusual Solution found\n");
			print_solution(solution);
		}
		return solution;
	}
	Solution get_random_solution()
	{
		Solution solution;
		int m, s;
		for (int i = 0; i < N; i++)
		{
			m = u(0, M);
			s = u(0, N);
			if (solution.value[m][s] == VALUE_DEFAULT)
				solution.value[m][s] = i;
			else
				i--;
		}
		if (!check_solution(solution))
		{
			printf("Warning: Unusual Solution found\n");
			print_solution(solution);
		}
		return solution;
	}
};


/*
2017-3-15
粒子群算法

Struct Pack:
编码存储结构

pack(): 将机器号、任务在机器内序列编码为int整数
unpack_machine_id(): 从包内解压出机器号
unpack_sequence_id(): 从包内解压出序列号

粒子群公式
v = r0 * v0 + r1 * (g - x) + r2 * (p - x)
x = x + v
*/

class PSOOptimizer
{
	// Strategy: assign jobs to machine, no order
	class Pack
	{
		int data[ALLOCATE_MAX_N];

	public:

		Pack() {}
		Pack(const Solution& solution)
		{
			initialize_from_solution(solution);
		}

		int pack_data(int machine_id, int order_in_sequence)
		{
			return machine_id * ALLOCATE_MAX_N + order_in_sequence;
		}

		int unpack_data_machine_id(int data) const
		{
			return data / ALLOCATE_MAX_N;
		}

		int unpack_data_order(int data) const
		{
			return data % ALLOCATE_MAX_N;
		}

		Solution get_solution() const
		{
			Solution result;

			for (int i = 0; i < N; i++)
				result.value[unpack_data_machine_id(data[i])][unpack_data_order(data[i])] = i;

			return result;
		}

		int get_fitness_value()
		{
			return evaluate(get_solution());
		}

		void initialize_as_sequence()
		{


			for (int i = 0; i < N; i++)
			{
				data[i] = pack_data(u(0, M), u(0, N));
				if (indexof(data, i, data[i]) != -1)
					i--;
			}



			//print_solution(get_solution());
			//printf("value = %d\n", get_fitness_value());

		}

		void initialize_as_binary()
		{
			for (int i = 0; i < N; i++)
			{
				data[i] = u(0, 100) < (SAVE_PROBABILITY * 100) ? 1 : 0;
			}
		}

		friend Pack operator+(const Pack& l, const Pack& r)
		{
			int index_1, index_2;
			Pack result;

			for (int i = 0; i < N; i++)
			{
				result.data[i] = l.data[i];
			}
			for (int i = 0; i < N; i++)
			{
				if (u(0, 100) < SELECTION_PROBABILITY * 100)
				{

					if (r.data[i] != VALUE_DEFAULT)
					{
						index_1 = indexof(result.data, N, r.data[i]);

						if (index_1 != -1)
						{
							index_2 = i;
							//if(index_1 == index_2) printf("...\n");
							swap(result.data[index_1], result.data[index_2]);
						}
						else
						{
							result.data[i] = r.data[i];
						}
					}
				}
			}
			return result;
		}

		friend Pack operator-(const Pack& l, const Pack& r)
		{
			Pack result;
			for (int i = 0; i < N; i++)
			{
				if (l.data[i] == r.data[i])
					result.data[i] = VALUE_DEFAULT;
				else
					result.data[i] = l.data[i];
			}
			return result;
		}

		friend Pack operator*(const Pack& l, const Pack& r)
		{
			Pack result;
			for (int i = 0; i < N; i++)
			{
				if (l.data[i] == 1)
					result.data[i] = r.data[i];
				else
					result.data[i] = VALUE_DEFAULT;
			}
			return result;
		}

		Pack& operator=(const Pack& r)
		{
			for (int i = 0; i < N; i++)
				data[i] = r.data[i];
			return *this;
		}

		void initialize_from_solution(const Solution& solution)
		{
			for (int i = 0; i < M; i++)
			{
				for (int j = 0; j < N; j++)
				{
					if (solution.value[i][j] != VALUE_DEFAULT)
					{
						data[solution.value[i][j]] = pack_data(i, j);
					}
				}
			}
		}

		void mutation()
		{

			for (int i = 0; i < N; i++)
			{
				if (u(0, 100) < MUTATION_PROBABILITY * 100)
				{
					int new_position = pack_data(u(0, M), u(0, N));
					int index = indexof(data, N, new_position);
					if (index == -1)
					{
						// new position is empty
						data[i] = new_position;
					}
					else
					{
						data[index] = data[i];
						data[i] = new_position;
						//printf("found\n");
					}
				}
			}

			//check_solution(get_solution());
		}
	};

	class Particle
	{
		Pack velocity;
		Pack position;

		Pack local_best;
		int local_best_value;

	public:

		void onPositionChanged(Pack& global_best, int& global_best_value)
		{
			int value = position.get_fitness_value();
			if (value < local_best_value)
			{
				local_best = position;
				local_best_value = value;

				if (value < global_best_value)
				{
					global_best = position;
					global_best_value = value;
				}
			}
		}

		void setPosition(Pack& new_position, Pack& global_best, int& global_best_value)
		{
			position = new_position;
			int value = position.get_fitness_value();

			onPositionChanged(global_best, global_best_value);
		}
		void update(Pack& global_best, int& global_best_value)
		{
			Pack r0, r1, r2;

			r0.initialize_as_binary();
			r1.initialize_as_binary();
			r2.initialize_as_binary();

			velocity = r0 * velocity + r1 * (local_best - position) + r2 * (global_best - position);
			velocity.mutation();

			position = position + velocity;

			onPositionChanged(global_best, global_best_value);

		}

		void initialize()
		{
			velocity.initialize_as_sequence();
			position.initialize_as_sequence();
			local_best.initialize_as_sequence();
			local_best_value = local_best.get_fitness_value();
		}
	};

	Pack global_best;
	Particle* particles;

	int global_best_value;
	int initial_spec_solution_count;

	void add_list_solution(const Solution& solution)
	{
		Pack pack(solution);
		particles[initial_spec_solution_count++].setPosition(pack, global_best, global_best_value);
	}

	void initialize()
	{
		for (int i = 0; i < POPULATION_SIZE; i++)
		{
			particles[i].initialize();
		}

		global_best.initialize_as_sequence();
		global_best_value = global_best.get_fitness_value();
		initial_spec_solution_count = 0;

		// add some list solution
		ListOptimizer listOptimizer;
		add_list_solution(listOptimizer.get_minc_solution());
		//printf("minc_added\n");
		add_list_solution(listOptimizer.get_miny_solution());
		//printf("miny_added\n");
		add_list_solution(listOptimizer.get_lpt_edd_solution());
		//printf("lpt_edd_added\n");
		add_list_solution(listOptimizer.get_spt_edd_solution());
		//printf("spt_edd_added\n");
		//add_list_solution(listOptimizer.get_min_makespan_first_solution());
		//add_list_solution(listOptimizer.get_min_proc_time_first_solution());

		tick_count = 0;
	}

	// 2017-5-3
	// 增加统计算法时间

	int tick_count;

public:

	PSOOptimizer()
	{
		particles = new Particle[POPULATION_SIZE];
	}

	~PSOOptimizer()
	{
		delete[] particles;
	}

	// 2017-5-3
	// 增加统计算法时间

	double get_running_time()
	{
		return tick_count / 1000.0;
	}

	Solution minimize()
	{
		initialize();

		int last_best_value = -1;
		int unchanged_iteration = 0, iteration = 0;

		// 2017-4-26
		// print data for graph
		char filename[256];
		sprintf(filename, "../Data/Pso%d-%d.txt", M, N);
		FILE *fp = NULL;
		if (PRINT_ITERATION_SOLUTION)
			fp = fopen(filename, "w");
		tick_count = GetTickCount();
		while (iteration < PSO_ITERATION_MAXN)
		{
			for (int i = 0; i < POPULATION_SIZE; i++)
			{
				particles[i].update(global_best, global_best_value);
			}
			if (PSO_UNCHANGED_ITERATION_MODE)
			{
				if (last_best_value != global_best_value)
				{
					last_best_value = global_best_value;
					//printf("new solution found, current value: %d\n", global_best_value);
					unchanged_iteration = 0;
				}
				else
				{
					unchanged_iteration++;
					if (unchanged_iteration > PSO_UNCHANGED_ITERATION_MAXN)
						break;
				}
			}
			iteration++;

			// 2017-4-26
			// print data for graph
			if (PRINT_ITERATION_SOLUTION)
				fprintf(fp, "%d\n", global_best_value);
		}

		tick_count = GetTickCount() - tick_count;

		if (PRINT_ITERATION_SOLUTION)
			fclose(fp);

		Solution current_solution = global_best.get_solution();

		if (!check_solution(current_solution))
		{
			printf("Warning: Unusual Solution found\n");
			print_solution(current_solution);
		}
		return current_solution;
	}

};

/*

2017-4-11
禁忌搜索算法

find_next_solution_random_method()
策略1:全排列搜索策略

find_next_solution_combination_method()
策略2:随机行走搜索策略

*/
class TabuOptimizer
{
	int recur;
	int tabu_pointer_1, tabu_pointer_2;

	// 2017-5-3
	// 增加统计算法时间
	int tick_count;

	bool terminated;
	// 2017-5-1
	// 添加短期记忆、中期记忆表
	Solution* tabu_list_1;
	Solution* tabu_list_2;

	bool solution_exists(const Solution& solution)
	{
		for (int i = 0; i < TABU_LIST_LENGTH; i++)
		{
			if (tabu_list_1[i].equals(solution))
				return true;
			if (tabu_list_2[i].equals(solution))
				return true;
		}
		return false;
	}

	void record_short_term_memory(const Solution& solution)
	{
		tabu_list_1[tabu_pointer_1++] = solution;
		if (tabu_pointer_1 >= TABU_LIST_LENGTH)
			tabu_pointer_1 = 0;
	}

	void record_medium_term_memory(const Solution& solution)
	{
		tabu_list_2[tabu_pointer_2++] = solution;
		if (tabu_pointer_2 >= TABU_LIST_LENGTH)
			tabu_pointer_2 = 0;
	}

	Solution find_neighbor_solution(const Solution& source)
	{
		Solution current = source;
		int m0 = u(0, M), m1 = u(0, M);
		int s0 = u(0, N), s1 = u(0, N);
		swap(current.value[m0][s0], current.value[m1][s1]);
		return current;
	}

	Solution find_next_solution_random_method(const Solution& source)
	{
		Solution current_solution = source;
		int current_solution_value = evaluate(current_solution), temp_value;
		bool updated = false;

		for (int i = 0; i < TABU_RANDOM_SELECTION_ITERATION; i++)
		{
			Solution next_solution = find_neighbor_solution(source);
			temp_value = evaluate(next_solution);

			// 2017-5-3
			// change search method

			if (!solution_exists(next_solution) && (temp_value < current_solution_value) || !updated)
			{
				record_short_term_memory(next_solution);
				current_solution_value = temp_value;
				current_solution = next_solution;
				updated = true;
			}
		}

		return current_solution;
	}

	Solution find_next_solution_combination_method(const Solution& source)
	{

		Solution local_best = source;
		Solution current = source;

		int current_value = evaluate(source);
		int local_best_value = current_value;
		bool updated = false;

		for (int i = 0; i < M; i++)
		{
			for (int i0 = 1; i0 < N; i0++)
			{
				for (int i1 = 0; i1 < i0; i1++)
				{
					if (source.value[i][i0] != VALUE_DEFAULT || source.value[i][i1] != VALUE_DEFAULT)
					{
						current = source;
						swap(current.value[i][i0], current.value[i][i1]);
						current_value = evaluate(current);

						if (current_value < local_best_value && !solution_exists(current))
						{
							local_best_value = current_value;
							local_best = current;
							updated = true;
						}
					}
				}
			}
		}

		for (int i = 0; i < N; i++)
		{
			for (int i0 = 1; i0 < M; i0++)
			{
				for (int i1 = 0; i1 < i0; i1++)
				{
					if (source.value[i0][i] != VALUE_DEFAULT || source.value[i1][i] != VALUE_DEFAULT)
					{
						current = source;
						swap(current.value[i0][i], current.value[i1][i]);
						current_value = evaluate(current);
						if (current_value < local_best_value && !solution_exists(current))
						{
							local_best_value = current_value;
							local_best = current;
							updated = true;
						}
					}
				}
			}
		}

		//terminated = !updated;

		return local_best;
	}

	void initialize()
	{
		tabu_pointer_1 = tabu_pointer_2 = 0;
		terminated = false;
		tick_count = 0;
	}

public:

	TabuOptimizer()
	{
		tabu_list_1 = new Solution[TABU_LIST_LENGTH];
		tabu_list_2 = new Solution[TABU_LIST_LENGTH];
	}

	~TabuOptimizer()
	{
		delete[] tabu_list_1;
		delete[] tabu_list_2;
	}

	Solution minimize()
	{
		initialize();
		ListOptimizer init_optimizer;


		Solution current_solution =
			init_optimizer.get_minc_solution();

		Solution best_solution = current_solution;

		tick_count = GetTickCount();
		int unchanged_iteration = 0, iteration = 0;
		int last_solution_value = 0, current_value = evaluate(current_solution), best_value = evaluate(best_solution);

		// 2017-4-26
		// print data for graph

		char filename[256];
		sprintf(filename, "../Data/Tabu%d-%d.txt", M, N);
		FILE *fp = NULL;
		if (PRINT_ITERATION_SOLUTION)
			fp = fopen(filename, "w");

		while (iteration < TABU_ITERATION_MAXN)
		{
			// 2017-5-3
			// change strategy of tabu search
			// record_solution(current_solution);
			current_solution = find_next_solution_random_method(current_solution);
			current_value = evaluate(current_solution);

			if (best_value > current_value)
			{
				record_medium_term_memory(current_solution);
				best_value = current_value;
				best_solution = current_solution;
			}

			if (TABU_UNCHANGED_ITERATION_MODE)
			{
				if (current_value != last_solution_value)
				{
					last_solution_value = current_value;
					unchanged_iteration = 0;
				}
				else
				{
					unchanged_iteration++;
					if (unchanged_iteration > TABU_UNCHANGED_ITERATION_MAXN)
						break;
				}
			}
			iteration++;

			if (terminated) break;

			// 2017-4-26
			// print data for graph
			if (PRINT_ITERATION_SOLUTION)
				fprintf(fp, "%d\n", current_value);
		}

		tick_count = GetTickCount() - tick_count;

		if (PRINT_ITERATION_SOLUTION)
			fclose(fp);

		if (!check_solution(current_solution))
		{
			printf("Warning: Unusual Solution found\n");
			print_solution(current_solution);
		}
		return best_solution;
	}
	// 2017-5-3
	// 增加统计算法时间
	double get_running_time()
	{
		return tick_count / 1000.0;
	}
};

void initialize()
{
	srand(time(0));
	if (GENERATE_JOB_ON_START)
		generate_test_data();
	load_test_data();
}

double calculate_improvement(int before_value, int after_value)
{
	// lower, better
	return 100.0 * (before_value - after_value) / before_value;
}

void create_record_file()
{
	FILE *fp = fopen("improvements_table.txt", "w");
	fprintf(fp, "M,N,b,PSO-LPT-EDD,PSO-SPT-EDD,PSO-MINC,PSO-MINY,PSO-T,TABU-LPT-EDD,TABU-SPT-EDD,TABU-MINC,TABU-MINY,TABU-T\n");
	fclose(fp);
}

// 2017-5-9
// 计算提升率
void append_record(double pso_value, double tabu_value, double pso_time, double tabu_time,
	double lpt_edd, double spt_edd, double minc, double miny)
{
	FILE *fp = fopen("improvements_table.txt", "a");
	// environment
	fprintf(fp, "%d,%d,%d,", M, N, b);
	// pso group
	double impr_pso_lpt_edd = calculate_improvement(lpt_edd, pso_value);
	double impr_pso_spt_edd = calculate_improvement(spt_edd, pso_value);
	double impr_pso_minc = calculate_improvement(minc, pso_value);
	double impr_pso_miny = calculate_improvement(miny, pso_value);
	fprintf(fp, "%.2f,%.2f,%.2f,%.2f,%.2f,", impr_pso_spt_edd, impr_pso_lpt_edd, impr_pso_minc, impr_pso_miny, pso_time);
	// tabu group
	double impr_tabu_lpt_edd = calculate_improvement(lpt_edd, tabu_value);
	double impr_tabu_spt_edd = calculate_improvement(spt_edd, tabu_value);
	double impr_tabu_minc = calculate_improvement(minc, tabu_value);
	double impr_tabu_miny = calculate_improvement(miny, tabu_value);
	fprintf(fp, "%.2f,%.2f,%.2f,%.2f,%.2f", impr_tabu_spt_edd, impr_tabu_lpt_edd, impr_tabu_minc, impr_tabu_miny, tabu_time);

	fprintf(fp, "\n");
	fclose(fp);
}
// 2017-4-26
// 引入批量测试

// 2017-5-3
// 增加统计算法时间
void batch_test()
{
	Solution solution;

	PSOOptimizer psoOptimizer;
	TabuOptimizer tabuOptimizer;
	ListOptimizer listOptimizer;

	double pso_value = 0, pso_value_sum = 0;
	double tabu_value = 0, tabu_value_sum = 0;

	double minc_value = 0, minc_value_sum = 0;
	double miny_value = 0, miny_value_sum = 0;
	double lpt_edd_value = 0, lpt_edd_value_sum = 0;
	double spt_edd_value = 0, spt_edd_value_sum = 0;

	double pso_time_value = 0, pso_time_value_sum = 0;
	double tabu_time_value = 0, tabu_time_value_sum = 0;

	FILE *fp = NULL;
	create_record_file();

	printf("M\tN\tb\tPSO\tPSO-t\tTABU\tTABU-t\tLPT-EDD\tSPT-EDD\tMINC\tMINY\n");
	if (PRINT_TABLE)
	{
		fp = fopen("output_table.txt", "w");
		fprintf(fp, "M\tN\tb\tPSO\tPSO-t\tTABU\tTABU-t\tLPT-EDD\tSPT-EDD\tMINC\tMINY\n");
	}

	for (int m_i = 0; m_i < sizeof(Ms) / sizeof(int); m_i++)
	{
		M = Ms[m_i];
		for (int n_i = 0; n_i < sizeof(NdMs) / sizeof(int); n_i++)
		{
			N = M * NdMs[n_i];
			for (int b_i = 0; b_i < sizeof(Bs) / sizeof(int); b_i++)
			{
				b = Bs[b_i];
				pso_time_value_sum = tabu_time_value_sum = pso_value_sum = tabu_value_sum = 0;
				lpt_edd_value_sum = spt_edd_value_sum = minc_value_sum = miny_value_sum = 0;
				for (int iteration = 0; iteration < ITERATION_IN_EACH_CASE; iteration++)
				{

					generate_test_data();
					printf("%d\t%d\t%d\t", M, N, b);
					// PSO
					solution = psoOptimizer.minimize();
					pso_value = evaluate(solution);
					pso_value_sum += pso_value;
					pso_time_value = psoOptimizer.get_running_time();
					pso_time_value_sum += pso_time_value;

					printf("%.0f\t", pso_value);
					printf("%.2f\t", pso_time_value);

					solution = tabuOptimizer.minimize();
					tabu_value = evaluate(solution);
					tabu_value_sum += tabu_value;
					tabu_time_value = tabuOptimizer.get_running_time();
					tabu_time_value_sum += tabu_time_value;

					printf("%.0f\t", tabu_value);
					printf("%.2f\t", tabu_time_value);

					solution = listOptimizer.get_lpt_edd_solution();
					lpt_edd_value = evaluate(solution);
					lpt_edd_value_sum += lpt_edd_value;

					printf("%.0f\t", lpt_edd_value);

					solution = listOptimizer.get_spt_edd_solution();
					spt_edd_value = evaluate(solution);
					spt_edd_value_sum += spt_edd_value;

					printf("%.0f\t", spt_edd_value);

					solution = listOptimizer.get_minc_solution();
					minc_value = evaluate(solution);
					minc_value_sum += minc_value;

					printf("%.0f\t", minc_value);

					solution = listOptimizer.get_miny_solution();
					miny_value = evaluate(solution);
					miny_value_sum += miny_value;

					printf("%.0f\n", miny_value);
				}

				pso_value = pso_value_sum / ITERATION_IN_EACH_CASE;
				tabu_value = tabu_value_sum / ITERATION_IN_EACH_CASE;
				pso_time_value = pso_time_value_sum / ITERATION_IN_EACH_CASE;
				tabu_time_value = tabu_time_value_sum / ITERATION_IN_EACH_CASE;

				lpt_edd_value = lpt_edd_value_sum / ITERATION_IN_EACH_CASE;
				spt_edd_value = spt_edd_value_sum / ITERATION_IN_EACH_CASE;
				minc_value = minc_value_sum / ITERATION_IN_EACH_CASE;
				miny_value = miny_value_sum / ITERATION_IN_EACH_CASE;

				if (PRINT_TABLE)
					fprintf(fp, "%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", M, N, b, pso_value, pso_time_value,
						tabu_value, tabu_time_value, lpt_edd_value, spt_edd_value, minc_value, miny_value);
				printf("%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t(AVG)\n", M, N, b, pso_value, pso_time_value,
					tabu_value, tabu_time_value, lpt_edd_value, spt_edd_value, minc_value, miny_value);
				append_record(pso_value, tabu_value, pso_time_value, tabu_time_value, lpt_edd_value, spt_edd_value, minc_value, miny_value);
			}
		}
	}
	if (PRINT_TABLE)
		fclose(fp);
}

int main()
{
	initialize();

	/*  2017-4-12 批量测试
	Solution solution;

	PSOOptimizer psoOptimizer;
	TabuOptimizer tabuOptimizer;
	ListOptimizer listOptimizer;

	solution = psoOptimizer.minimize();
	print_solution(solution, "Particle Swarm Optimization");

	solution = tabuOptimizer.minimize();
	print_solution(solution, "Tabu Search");

	solution = listOptimizer.get_min_makespan_first_solution();
	print_solution(solution, "List: machineID = argmin(makespan)");

	solution = listOptimizer.get_min_proc_time_first_solution();
	print_solution(solution, "List: machineID = argmin(process time on machine)");
	*/
	batch_test();

	//printf("list result: %d\n", evaluate(solution));
	return 0;
}
