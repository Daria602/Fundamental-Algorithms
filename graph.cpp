#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <algorithm>
#include <stack>
#include <tuple>

using namespace std;
# define INF 0x3f3f3f3f

ifstream fin("test.txt");

deque<int> count_sort_reverse(deque<int> to_sort)
{
	vector<int> freq_list;
	deque<int> result;
	freq_list.resize(*max_element(to_sort.begin(), to_sort.end()) + 1);
	for (unsigned int i = 0; i < to_sort.size(); i++)
	{
		freq_list[to_sort[i]]++;
	}
	for (unsigned int i = 0; i < freq_list.size(); i++)
	{
		while (freq_list[i] > 0)
		{
			result.push_front(i);
			freq_list[i]--;
		}
	}
	return result;
}

//Graph class, all the nodes start from 0, modify the add_edge in main() and all the "cout"s in member functions to get from 1
class Graph
{
private:
	//number of nodes
	int n;
	//number of edges
	int e;
	//bool if graph is oriented
	bool oriented;
	//adj list for graph representation
	vector<vector<int>> adj_list;
public:
	Graph() {}
	Graph(int n, bool oriented) 
	{
		this->n = n;
		this->e = 0;
		this->oriented = oriented;
		//populate adj list with empty vectors
		vector<int> empty;
		for (int i = 0; i < n; i++)
		{
			this->adj_list.push_back(empty);
		}
	}
	virtual ~Graph() {}

	void add_edge(int node1, int node2)
	{
		this->e++;
		this->adj_list[node1].push_back(node2);

		//if graph is not oriented, then push from the second node
		if (!this->oriented)
		{
			this->adj_list[node2].push_back(node1);
		}
	}

private:
	//dfs traversal
	void dfs(int start_node, vector<bool>& visited)
	{
		visited[start_node] = true;
		cout << start_node << " ";

		for (unsigned int i = 0; i < adj_list[start_node].size(); i++)
		{
			int current_node;
			current_node = adj_list[start_node][i];
			if (visited[current_node] == false)
			{
				dfs(current_node, visited);
			}
		}

	}

	void find_biconnected_comp(int current_node, int parent_node, int current_distance,
		vector<int>& distance, vector<int>& shortest_distance, vector<bool>& visited, vector<int>& nodes, vector<vector<int>>& bicon_comps)
	{
		distance[current_node] = current_distance;
		shortest_distance[current_node] = current_distance;
		visited[current_node] = true;
		nodes.push_back(current_node);

		for (unsigned int i = 0; i < this->adj_list[current_node].size(); i++)
		{
			int next_node;
			next_node = this->adj_list[current_node][i];
			if (next_node != parent_node)
			{
				if (!visited[next_node])
				{
					find_biconnected_comp(next_node, current_node, current_distance + 1, distance, shortest_distance, visited, nodes, bicon_comps);

					shortest_distance[current_node] = min(shortest_distance[current_node], shortest_distance[next_node]);

					if (distance[current_node] <= shortest_distance[next_node])
					{

						vector<int> bicon_comp;
						bicon_comp.push_back(current_node);
						int node = nodes.back();
						nodes.pop_back();
						bicon_comp.push_back(node);

						while (node != next_node)
						{
							node = nodes.back();
							nodes.pop_back();
							bicon_comp.push_back(node);

						}

						bicon_comps.push_back(bicon_comp);

					}
				}
				else
				{

					//modify the shortest distance if the node is already visited
					shortest_distance[current_node] = min(shortest_distance[current_node], distance[next_node]);
				}
			}

		}
	}

	void sccdfs(int current_node, int& id, stack<int>& stck, vector<bool>& on_stack, vector<int>& ids, vector<int>& low_link_values, vector<vector<int>>& str_con_comps)
	{
		stck.push(current_node);
		on_stack[current_node] = true;
		ids[current_node] = id;
		low_link_values[current_node] = id;
		id++;

		for (unsigned int i = 0; i < this->adj_list[current_node].size(); i++)
		{
			int neighbour;
			neighbour = adj_list[current_node][i];

			if (ids[neighbour] == -1)
			{
				sccdfs(neighbour, id, stck, on_stack, ids, low_link_values, str_con_comps);
			}
			if (on_stack[neighbour])
			{
				low_link_values[current_node] = min(low_link_values[current_node], low_link_values[neighbour]);
			}
		}

		if (ids[current_node] == low_link_values[current_node])
		{
			vector<int> str_con_comp;
			int node = stck.top();
			while (node != current_node)
			{
				str_con_comp.push_back(node);
				stck.pop();
				on_stack[node] = false;
				low_link_values[node] = ids[current_node];
				node = stck.top();
			}
			str_con_comp.push_back(current_node);
			stck.pop();

			str_con_comps.push_back(str_con_comp);
		}

	}

	void top_sort_dfs(int current_node, vector<bool>& visited, stack<int>& stck)
	{
		visited[current_node] = true;

		for (int i = this->adj_list[current_node].size() - 1; i >= 0; i--)
		{
			int neighbour;
			neighbour = adj_list[current_node][i];
			if (!visited[neighbour])
			{
				top_sort_dfs(neighbour, visited, stck);
			}

		}
		stck.push(current_node);
	}

	void ccdfs(int current_node, int parent_node, int& distance, vector<bool>& visited, vector<int>& distances, vector<vector<int>>& result)
	{
		visited[current_node] = true;
		distances[current_node] = distance++;
		int current_distance = distances[current_node];

		for (int neighbour : adj_list[current_node])
		{
			if (neighbour == parent_node) continue;
			if (!visited[neighbour])
			{
				ccdfs(neighbour, current_node, distance, visited, distances, result);
			}

			distances[current_node] = min(distances[current_node], distances[neighbour]);

			if (current_distance < distances[neighbour])
			{
				vector<int> pair;
				pair.push_back(current_node);
				pair.push_back(neighbour);

				result.push_back(pair);
			}
		}

	}

	pair<int, int> bfs_darb(int start_node)
	{
		deque<int> unvisited;
		vector<int> visited;
		vector<int> answer;
		vector<bool> unvisited2(this->n, false);
		vector<bool> visited2(this->n, false);

		for (int i = 0; i < this->n; i++)
		{
			answer.push_back(-1);
		}

		unvisited.push_back(start_node);
		unvisited2[start_node] = true;
		answer[start_node] = 1;
		int last_node = 0;

		while (!unvisited.empty())
		{
			for (unsigned int i = 0; i < this->adj_list[start_node].size(); i++)
			{
				int key;
				key = adj_list[start_node][i];

				if (!visited2[key] && !unvisited2[key])
				{
					unvisited.push_back(key);
					unvisited2[key] = true;
					answer[key] = answer[start_node] + 1;
					last_node = key;

				}

			}

			visited.push_back(start_node);
			visited2[start_node] = true;
			unvisited.pop_front();
			if (!unvisited.empty())
			{
				start_node = unvisited[0];
			}
		}
		int diameter = answer[last_node];
		return make_pair(last_node, diameter);
	}


public:

	//general bfs traversal (also works for graph that is not connected)
	void start_bfs(int start_node)
	{
		vector<int> visited;
		for (int i = 0; i < this->n; i++)
		{
			if (find(visited.begin(), visited.end(), i) == visited.end())
			{
				bfs(i, visited);
			}
		}

		for (unsigned int i = 0; i < visited.size(); i++)
		{
			cout << visited[i] << " ";
		}
	}
	void bfs(int start_node, vector<int> &visited)
	{
		deque<int> unvisited;
		unvisited.push_back(start_node);
		while (!unvisited.empty())
		{
			for (unsigned int i = 0; i < this->adj_list[start_node].size(); i++)
			{
				int key;
				key = adj_list[start_node][i];

				bool inVisited, inUnvisited;
				inVisited = find(visited.begin(), visited.end(), key) != visited.end();
				inUnvisited = find(unvisited.begin(), unvisited.end(), key) != unvisited.end();

				if (!inVisited && !inUnvisited)
				{
					unvisited.push_back(key);
				}
			}

			visited.push_back(start_node);
			unvisited.pop_front();
			if (!unvisited.empty())
			{
				start_node = unvisited[0];
			}

		}

	}
	//bfs for infoarena problem1
	// test (score 100) https://www.infoarena.ro/job_detail/2792989
	void bfsinfoarena(int start_node)
	{
		deque<int> unvisited;
		vector<int> visited;
		vector<int> answer;
		vector<bool> unvisited2(this->n, false);
		vector<bool> visited2(this->n, false);

		for (int i = 0; i < this->n; i++)
		{
			answer.push_back(-1);
		}

		unvisited.push_back(start_node);
		unvisited2[start_node] = true;
		answer[start_node] = 0;


		while (!unvisited.empty())
		{
			for (unsigned int i = 0; i < this->adj_list[start_node].size(); i++)
			{
				int key;
				key = adj_list[start_node][i];

				if (!visited2[key] && !unvisited2[key])
				{
					unvisited.push_back(key);
					unvisited2[key] = true;
					answer[key] = answer[start_node]+1;
				}

			}

			visited.push_back(start_node);
			visited2[start_node] = true;
			unvisited.pop_front();
			if (!unvisited.empty())
			{
				start_node = unvisited[0];
			}
		}

		for (unsigned int i = 0; i < answer.size(); i++)
		{
			cout << answer[i] << " ";
		}
		
	}
	
	//general dfs traversal 
	void start_dfs(int start_node)
	{
		vector<bool> visited;
		for (int i = 0; i < this->n; i++)
		{
			visited.push_back(false);
		}

		dfs(start_node, visited);

		//if graph is not connected, this will show the rest
		for (int i = 0; i < this->n; i++)
		{
			if (!visited[i])
			{
				dfs(i, visited);
			}
		}
	}


	//to see the number of connected elements in the graph for infoarena problem2
	//test (score 100) https://www.infoarena.ro/job_detail/2792288 
	int start_dfsinfoarena()
	{
		vector<bool> visited;
		int cnt = 0;
		for (int i = 0; i < this->n; i++)
		{
			visited.push_back(false);
		}

		for (int i = 0; i < this->n; i++)
		{
			if (!visited[i])
			{
				cnt++;
				dfs(i, visited);
			}
		}
		return cnt;
	}

	//biconnected components problem3
	//test (score 90, time limit) https://infoarena.ro/job_detail/2792556
	tuple<int, vector<vector<int>>> biconnected_components()
	{
		//setting up all the vectors
		vector<vector<int>> bicon_comps;
		vector<int> nodes;
		vector<bool> visited;
		vector<int> distance;
		vector<int> shortest_distance;

		for (int i = 0; i < this->n; i++)
		{
			visited.push_back(false);
			distance.push_back(-1);
			shortest_distance.push_back(-1);
		}

		find_biconnected_comp(0, -1, 0, distance, shortest_distance, visited, nodes, bicon_comps);

		int nr_comp_bicon = bicon_comps.size();
		tuple<int, vector<vector<int>>> result;
		result = make_tuple(nr_comp_bicon, bicon_comps);
		return result;
	}


	//strongly connected components problem4
	//test (score 90, time limit) https://infoarena.ro/job_detail/2792939
	tuple<int, vector<vector<int>>> strongly_connected_components()
	{
		//used to give each node an id
		int id = 0;
		int scc_count;

		vector<int> ids(this->n, -1);
		vector<int> low_link_values(this->n, 0);
		vector<bool> on_stack(this->n, false);
		stack<int> stck;
		vector<vector<int>> str_con_comps;

		for (int i = 0; i < this->n; i++)
		{
			if (ids[i] == -1)
			{
				sccdfs(i, id, stck, on_stack, ids, low_link_values, str_con_comps);
			}
		}

		if (!stck.empty())
		{
			scc_count = stck.size() + str_con_comps.size();
		}
		else
		{
			scc_count = str_con_comps.size();
		}
		

		while (!stck.empty())
		{
			vector<int> singular_comp;
			int nd = stck.top();
			singular_comp.push_back(nd);
			str_con_comps.push_back(singular_comp);
			stck.pop();
		}

		tuple<int, vector<vector<int>>> result;
		result = make_tuple(scc_count, str_con_comps);

		return result;


	}


	//topological sort problem5
	//test (score 100) https://infoarena.ro/job_detail/2792978
	vector<int> topological_sort()
	{
		vector<bool> visited(this->n, false);
		stack<int> stck;

		for (int i = 0; i < this->n; i++)
		{
			if (!visited[i])
			{
				top_sort_dfs(i, visited, stck);
			}
		}
		vector<int> result;

		while (!stck.empty())
		{
			result.push_back(stck.top());
			stck.pop();
		}

		return result;
	}
	
	// find if a degree sequence can form a simple graph problem6
	bool havel_hakimi()
	{
		int N;
		fin >> N;
		deque<int> degrees;
		for (int i = 0; i < N; i++)
		{
			int degree;
			fin >> degree;
			degrees.push_back(degree);
		}

		while (true)
		{
			degrees = count_sort_reverse(degrees);

			if (degrees[0] == 0)
			{
				return true;
			}

			int first_element = degrees.front();
			degrees.pop_front();

			if (first_element > degrees.size())
			{
				return false;
			}

			for (int i = 0; i < first_element; i++)
			{
				degrees[i]--;

				if (degrees[i] < 0)
				{
					return false;
				}
			}
		}
	}
	
	//critical connections problem7
	//test leetcode Runtime: 608 ms
	vector<vector<int>> critical_connections()
	{
		vector<int> empty;
		vector<vector<int>> result;


		int distance = 0;
		vector<bool> visited(this->n, false);
		vector<int> distances(this->n, -1);

		ccdfs(0, -1, distance, visited, distances, result);

		return result;
	}

	//minimum spanning tree algorithm 
	// test (score 100) https://infoarena.ro/job_detail/2802275
	tuple<int, int, vector<tuple<int, int>>> apm()
	{
		vector<vector<tuple<int, int>>> adj_list;
		vector<tuple<int, int>> empty;
		vector<bool> visited;
		priority_queue<tuple<int, int, int>, vector<tuple<int, int, int>>, greater<tuple<int, int, int>> > pq;
		vector<tuple<int, int>> answer;
		int N, E;
		fin >> N >> E;

		for (int i = 0; i < N; i++)
		{
			adj_list.push_back(empty);
			visited.push_back(false);
		}

		int edge1, edge2, cost;
		for (int i = 0; i < E; i++)
		{
			fin >> edge1 >> edge2 >> cost;
			edge1--;
			edge2--;
			adj_list[edge1].push_back(make_tuple(edge2, cost));
			adj_list[edge2].push_back(make_tuple(edge1, cost));
		}

		int sum_apm = 0;
		int optimized_nr_edges_answer = N - 1;
		int optimized_nr_edges = N - 1;


		int start_node = 0;

		for (unsigned int i = 0; i < adj_list[start_node].size(); i++)
		{
			pq.push(make_tuple(get<1>(adj_list[start_node][i]), start_node, get<0>(adj_list[start_node][i])));
		}

		visited[start_node] = true;

		while (!pq.empty() && optimized_nr_edges > 0)
		{
			tuple <int, int, int> mn;
			mn = pq.top();
			while (visited[get<2>(mn)])
			{
				pq.pop();
				mn = pq.top();
			}
			sum_apm = sum_apm + get<0>(mn);
			answer.push_back(make_tuple(get<1>(mn), get<2>(mn)));
			pq.pop();
			optimized_nr_edges--;
			start_node = get<2>(mn);

			for (unsigned int i = 0; i < adj_list[start_node].size(); i++)
			{
				if (!visited[get<0>(adj_list[start_node][i])])
				{
					pq.push(make_tuple(get<1>(adj_list[start_node][i]), start_node, get<0>(adj_list[start_node][i])));
				}

			}
			visited[start_node] = true;
		}
		tuple<int, int, vector<tuple<int, int>>> result;
		result = make_tuple(sum_apm, optimized_nr_edges_answer, answer);

		return result;


	}


	//dijkstra's algorithm to find the shortest path from one vertex to the others
	//test (score 90) https://infoarena.ro/job_detail/2813068
	vector<int> dijkstra()
	{
		vector<vector<tuple<int, int>>> adj_list;
		vector<tuple<int, int>> empty;
		vector<bool> visited;
		vector<int> distances;

		//reading the values
		int N, E;
		fin >> N >> E;

		for (int i = 0; i < N; i++)
		{
			adj_list.push_back(empty);
			visited.push_back(false);
			distances.push_back(INF);
		}

		int edge1, edge2, cost;
		for (int i = 0; i < E; i++)
		{
			fin >> edge1 >> edge2 >> cost;
			edge1--;
			edge2--;
			adj_list[edge1].push_back(make_tuple(edge2, cost));
			//adj_list[edge2].push_back(make_tuple(edge1, cost));
		}

		//the algorithm
		//min-heap sort by cost
		priority_queue<tuple<int, int>, vector<tuple<int, int>>, greater<tuple<int, int>> > pq;
		//distance from the first node to itself is 0
		// push initial node to heap
		int initial_node = 0;
		distances[initial_node] = 0;
		pq.push(make_tuple(0, initial_node));

		//pop min from heap after visiting all the neighbours
		while (!pq.empty())
		{
			int node;
			node = get<1>(pq.top());
			pq.pop();
			visited[node] = true;

			for (tuple<int, int> neighbour : adj_list[node])
			{
				if (visited[get<0>(neighbour)])
				{
					continue;
				}

				int new_dist = distances[node] + get<1>(neighbour);

				if (new_dist < distances[get<0>(neighbour)])
				{
					distances[get<0>(neighbour)] = new_dist;
					if (!visited[get<0>(neighbour)])
					{
						pq.push(make_tuple(new_dist, get<0>(neighbour)));
					}
				}
			}

		}


		for (unsigned int i = 1; i < distances.size(); i++)
		{
			if (distances[i] == INF)
			{
				distances[i] = 0;
			}
		}

		return distances;
	}


	//floyd-warshall algorithm to find min distance from every node to the other
	//test (score 100) https://infoarena.ro/job_detail/2812967
	vector<vector<int>> royfloyd()
	{
		vector<vector<int>> distanceMatrix(101);

		int N;
		fin >> N;
		this->n = N;
		for (int i = 0; i < this->n; i++)
		{
			for (int j = 0; j < this->n; j++)
			{
				int distance;
				fin >> distance;
				if (i != j && distance == 0)
				{
					distance = INF;
				}
				distanceMatrix[i].push_back(distance);
			}
		}



		for (int k = 0; k < this->n; k++)
		{
			for (int i = 0; i < this->n; i++)
			{
				for (int j = 0; j < this->n; j++)
				{
					distanceMatrix[i][j] = min(distanceMatrix[i][j], distanceMatrix[i][k] + distanceMatrix[k][j]);
				}
			}
		}
		vector<vector<int>> result;
		vector<int> empty;
		for (int i = 0; i < this->n; i++)
		{
			result.push_back(empty);
		}

		for (int i = 0; i < this->n; i++)
		{
			for (int j = 0; j < this->n; j++)
			{
				if (distanceMatrix[i][j] == INF)
				{
					distanceMatrix[i][j] = 0;
				}
				result[i].push_back(distanceMatrix[i][j]);
			}
		}

		return result;
	}


	// bfs traversal to find the diameter of n-ary tree (graph)
	// test (score 100) https://infoarena.ro/job_detail/2813142
	int darb()
	{


		//starting node
		int SN = 1;

		int n1, n2;
		for (int i = 0; i < this->n; i++)
		{
			fin >> n1 >> n2;
			add_edge(n1 - 1, n2 - 1);
		}

		int last_el = bfs_darb(SN - 1).first;
		return bfs_darb(last_el).second;

	}
	// find eulerian cycle
	// test (score 100) https://infoarena.ro/job_detail/2820205
	vector<int> euler_cycle()
	{
		fin >> this->n;

		//created adj list with all the edges being indexed
		vector<pair<int, int>> empty;
		vector<vector<pair<int, int>>> adj_list_indexed;
		for (int i = 0; i < this->n; i++)
		{
			adj_list_indexed.push_back(empty);
		}

		fin >> this->e;
		for (int i = 0; i < this->e; i++)
		{
			int n1, n2;
			fin >> n1;
			n1--;
			fin >> n2;
			n2--;

			adj_list_indexed[n1].push_back(make_pair(n2, i));
			adj_list_indexed[n2].push_back(make_pair(n1, i));

		}

		vector<int> res;

		for (int i = 0; i < adj_list_indexed.size(); i++)
		{
			if (adj_list_indexed[i].size() % 2 != 0)
			{

				res.push_back(-1);
				res.push_back(-1);
				return res;
			}
		}

		//bool vector to see if edge was already removed
		vector<bool> edge_removed(this->e, false);

		//create stack to remember the path
		stack<int> path;
		//for cycle we can start from any node, so let's start with 0
		path.push(0);

		while (!path.empty())
		{
			int current_node = path.top();

			if (!adj_list_indexed[current_node].empty())
			{
				pair<int, int> next_node_index = adj_list_indexed[current_node].back();
				adj_list_indexed[current_node].pop_back();

				if (!edge_removed[next_node_index.second])
				{
					edge_removed[next_node_index.second] = true;
					path.push(next_node_index.first);
				}
			}
			else
			{
				res.push_back(current_node + 1);
				path.pop();
			}
		}

		return res;

	}
	//find hamilton cycle with min cost
	//test (score 100) https://infoarena.ro/job_detail/2820094
	int hamilton_cycle_min_cost()
	{
		fin >> this->n;
		vector<vector<pair<int, int> >> adj_list_with_cost;
		vector<pair<int, int>> empty;
		for (int i = 0; i < this->n; i++)
		{
			adj_list_with_cost.push_back(empty);
		}
		int edges;
		fin >> edges;
		this->e = edges;
		for (int i = 0; i < this->e; i++)
		{
			int n1, n2, cost;
			fin >> n1 >> n2 >> cost;
			adj_list_with_cost[n1].push_back(make_pair(n2, cost));

		}

		int answer = INF;
		int nr_nodes = 1 << this->n;


		vector<int> v(this->n, INF);
		vector<vector<int>> costs(nr_nodes, v);

		costs[1][0] = 0;

		for (int i = 0; i < nr_nodes; i++)
		{
			for (int j = 0; j < this->n; j++)
			{
				if (i & (1 << j))
				{
					for (int k = 0; k < adj_list_with_cost[j].size(); k++)
					{
						if (i & (1 << adj_list_with_cost[j][k].first))
						{
							costs[i][j] = min(costs[i][j], costs[i ^ (1 << j)][adj_list_with_cost[j][k].first] + adj_list_with_cost[j][k].second);
						}
					}
				}
			}
		}

		for (int i = 0; i < adj_list_with_cost[0].size(); i++)
		{
			answer = min(answer, costs[nr_nodes - 1][adj_list_with_cost[0][i].first] + adj_list_with_cost[0][i].second);
		}
		return answer;


	}

};

pair<Graph,int> read_with_starting_node(bool orient)
{
	//number of nodes
	int N;
	//number of edges
	int E;
	//starting node
	int SN;
	fin >> N >> E;
	
	fin >> SN;

	Graph graph(N, orient);
	int n1, n2;
	for (int i = 0; i < E; i++)
	{
		fin >> n1 >> n2;
		graph.add_edge(n1, n2);
	}

	return make_pair(graph,SN);
}

Graph read_just_the_connections(bool orient)
{
	//number of nodes
	int N;
	//number of edges
	int E;
	//starting node
	int SN = 0;
	fin >> N >> E;
	//for problems that give starting node
	//fin >> SN

	//true if graph is oriented
	Graph graph(N, true);
	int n1, n2;
	for (int i = 0; i < E; i++)
	{
		fin >> n1 >> n2;
		graph.add_edge(n1, n2);
	}
	return graph;
}


int main()
{
	//bool for graph, true if oriented
	bool orient = false;

	//for bfsinfoarena 
	/*
	pair<Graph, int>  graph_and_start_node;
	graph_and_start_node = read_with_starting_node(orient);

	Graph graph = graph_and_start_node.first;
	graph.bfsinfoarena(graph_and_start_node.second);
	*/
	
	//for the problems below
	//Graph graph;
	//graph = read_just_the_connections(orient);

	//int result_nr_comp_conexe = graph.start_dfsinfoarena();
	//tuple<int, vector<vector<int>>> result_bicon_comp = graph.biconnected_components();
	//tuple<int, vector<vector<int>>> result_scc = graph.strongly_connected_components();
	//vector<int> result_top_sort = graph.topological_sort();
	//vector<vector<int>> result_cc = graph.critical_connections();

	//for the rest

	
	//Graph graph;
	//bool result_havel_hakimi = graph.havel_hakimi();
	//tuple<int, int, vector<tuple<int, int>>> result_apm = graph.apm();
	
	
	//vector<int> result_dijkstra = graph.dijkstra();

	//vector<vector<int>> result_royfloyd = graph.royfloyd();

	//vector<int> result_euler_cycle = graph.euler_cycle();
	//int result_ham = graph.hamilton_cycle_min_cost();

	//for diameter of n-ary tree
	/*
	int N;
	fin >> N;
	Graph g(N, false);
	int result_darb = g.darb();
	*/

	return 0;
	
}
