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
	
	

	//dfs traversal
	void dfs(int start_node, vector<bool> &visited)
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
	void start_dfsinfoarena()
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
		cout << endl;
		cout << "Nr de componente conexe: " << cnt;
	}

	//biconnected components problem3
	//test (score 90, time limit) https://infoarena.ro/job_detail/2792556

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

	void biconnected_components()
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

		cout << "Nr of biconnected components: " << bicon_comps.size() << endl;

		for (unsigned int i = 0; i < bicon_comps.size(); i++)
		{
			for (unsigned int k = 0; k < bicon_comps[i].size(); k++)
			{
				cout << bicon_comps[i][k] << " ";
			}

			cout << endl;
		}
	}
	//strongly connected components problem4
	//test (score 90, time limit) https://infoarena.ro/job_detail/2792939
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

	void strongly_connected_components()
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

		cout << scc_count << "\n";
		for (unsigned int i = 0; i < str_con_comps.size(); i++)
		{
			for (unsigned int k = 0; k < str_con_comps[i].size(); k++)
			{
				cout << str_con_comps[i][k]<< " ";
			}
			cout << "\n";
		}
		while (!stck.empty())
		{
			int nd = stck.top();
			cout << nd << "\n";
			stck.pop();
		}


	}
	//topological sort problem5
	//test (score 100) https://infoarena.ro/job_detail/2792978
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

	void topological_sort()
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

		while (!stck.empty())
		{
			cout << stck.top() << " ";
			stck.pop();
		}
	}
	// find if a degree sequence can form a simple graph problem6
	bool havel_hakimi()
	{
		int N;
		fin >> N;
		Graph graph(N, true); // just to use the functions from the class
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
				cout << "true";
				return true;
			}

			int first_element = degrees.front();
			degrees.pop_front();

			if (first_element > degrees.size())
			{
				cout << "false";
				return false;
			}

			for (int i = 0; i < first_element; i++)
			{
				degrees[i]--;

				if (degrees[i] < 0)
				{
					cout << "false";
					return false;
				}
			}
		}
	}
	
	//critical connections problem7
	//test leetcode Runtime: 608 ms
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
	void critical_connections()
	{
		vector<int> empty;
		vector<vector<int>> result;


		int distance = 0;
		vector<bool> visited(this->n, false);
		vector<int> distances(this->n, -1);

		ccdfs(0, -1, distance, visited, distances, result);

		for (unsigned int i = 0; i < result.size(); i++)
		{
			cout << result[i][0] << " " << result[i][1] << endl;
		}
	}

	//minimum spanning tree algorithm 
	// test (score 100) https://infoarena.ro/job_detail/2802275
	void apm()
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
		cout << sum_apm << "\n" << optimized_nr_edges_answer << "\n";

		for (unsigned int i = 0; i < answer.size(); i++)
		{
			cout << get<0>(answer[i]) << " " << get<1>(answer[i])  << "\n";
		}


	}
	//dijkstra's algorithm to find the shortest path from one vertex to the others
	//test (score 90) https://infoarena.ro/job_detail/2802388
	void dijkstra()
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
		priority_queue<tuple<int, int>, vector<tuple<int, int>>, greater<tuple<int, int>> > pq;
		distances[0] = 0;
		pq.push(make_tuple(0, 0));

		while (!pq.empty())
		{
			int node, dist;
			node = get<1>(pq.top());
			dist = get<1>(pq.top());
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
					pq.push(make_tuple(new_dist, get<0>(neighbour)));
				}
			}

		}


		for (unsigned int i = 1; i < distances.size(); i++)
		{
			if (distances[i] == INF)
			{
				cout << "0 ";
			}
			else
			{
				cout << distances[i] << " ";
			}
		}
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

	//graph.start_dfsinfoarena();
	//graph.biconnected_components();
	//graph.strongly_connected_components();
	//graph.topological_sort();
	//graph.critical_connections();

	//for havel-hakimi, minimum spanning tree (apm), dijkstra

	
	Graph graph;
	//graph.havel_hakimi();
	//graph.apm();
	//graph.dijkstra();

	return 0;
	
}
