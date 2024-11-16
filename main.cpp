#include <algorithm>
#include <iostream>
#include <random>
#include <fstream>
#include <set>
#include <unordered_set>
#include <vector>

using namespace std;

struct Statistics {
	vector<double> data;
	double mean, dev, total;
	Statistics(vector<double>&& data_): data(data_), dev(0) {
		sort(data.begin(), data.end());
		total=double(reduce(data.begin(), data.end()));
		mean=total/double(data.size());

		for (double x: data) dev+=(x-mean)*(x-mean);
		dev/=data.size();
		dev = sqrt(dev);
	}

	friend ostream& operator<<(ostream& os, Statistics const& bench) {
		cout<<bench.data.size()<<" values, "<<bench.mean<<" mean, "<<bench.dev<<" stddev"<<endl;
		for (size_t i=0; i<=4; i++) {
			cout<<"Q"<<i<<": "<<bench.data[(i*(bench.data.size()-1))/4];
			i==4 ? cout<<endl : cout<<", ";
		}

		return os;
	}

	void plot(ostream& os) {
		os<<"set boxwidth 0.5"<<endl;
		os<<"set style fill solid"<<endl;
		os<<"set term pdfcairo"<<endl;
		os<<"set output \"./bench.pdf\""<<endl;
		os<<"set yrange [0:]"<<endl;
		os<<"set xrange [0:]"<<endl;
		os<<"plot '-' using 1 bins notitle with boxes"<<endl;

		for (size_t i=0; i<data.size(); i++) {
			os<<data[i]<<endl;
		}

		os<<"e"<<endl;
	}
};

int main(int argc, char** argv) {
	bool stay_home=false;
	if (argc>1 && string(argv[1])=="--stay-home") {
		stay_home=true;
		cout<<"accounting for staying home while sick\n";
	}

	constexpr int n = 2000, m = 20;
	constexpr double rewire_prob = 0.15, initial_infection_prob=0.1;
	constexpr int K=10;

	constexpr int num_days = 100;
	constexpr int time_steps_in_day = 50;
	constexpr int plot_every_time_step = 50;
	constexpr double day_person_infection_rate = 0.03;
	constexpr double day_person_recovery_rate = 0.2;
	constexpr int num_trials = 10;

	vector<int> total_infected(num_days*time_steps_in_day/plot_every_time_step, 0);
	vector<int> total_susceptible=total_infected;
	for (int trial=1; trial<=num_trials; trial++) {
		mt19937 rng(1234+trial);

		vector<unordered_set<int>> edges(n);
		for (int i=0; i<n; i++) {
			for (int j=1; j<=m/2; j++) {
				int edge_to=-1;
				if (uniform_real_distribution<>()(rng)<rewire_prob) {
					edge_to = (uniform_int_distribution<>(0,n-m-2)(rng) + i+m/2+1)%n;
					if (edges[i].contains(edge_to) || (edge_to<i && edges[edge_to].contains(i)))
						edge_to=-1;
				}

				if (edge_to==-1) {
					edge_to = j+i;
					if (edge_to>=n) edge_to-=n;
				}

				edges[i].insert(edge_to);
				edges[edge_to].insert(i);
			}
		}

		cout<<"Computing graph statistics ("<<n*m<<" edges)\n";

		vector<double> clustering_coefficients, degree, path_length;
		vector<vector<int>> apsp(n, vector<int>(n, INT_MAX));
		for (int i=0; i<n; i++) {
			degree.push_back(edges[i].size());

			int closed_tris=0, num_triples=edges[i].size()*(edges[i].size()-1)/2;
			for (int j: edges[i]) {
				apsp[i][j]=1;
				for (int k: edges[i]) {
					if (j>k && edges[j].contains(k)) closed_tris++;
				}
			}

			clustering_coefficients.push_back(1.0*closed_tris/num_triples);
		}

		for (int i=0; i<n; i++) {
			queue<array<int,2>> bfs;
			vector<bool> visit(n);
			bfs.push({i,0});
			visit[i]=true;
			while (bfs.size()) {
				auto [x,l] = bfs.front();
				bfs.pop();
				for (int y: edges[x]) {
					if (!visit[y]) {
						if (y<i) path_length.push_back(l+1);
						bfs.push({y, l+1});
						visit[y]=true;
					}
				}
			}
		}

		Statistics cluster_stat(std::move(clustering_coefficients));
		Statistics path_stat(std::move(path_length));
		cout<<"Clustering coefficients:\n"<<cluster_stat<<"\n";
		cout<<"Vertex degree:\n"<<Statistics(std::move(degree))<<"\n";
		cout<<"Path length:\n"<<path_stat<<"\n";
		cout<<"Small world coefficient: "<< cluster_stat.mean/(1.0*m/n) / (path_stat.mean/(log(n)/log(m))) <<"\n";

		double time_step_infection_rate = 1-pow(1-day_person_infection_rate, 1.0/time_steps_in_day);
		double time_step_recovery_rate = 1-pow(1-day_person_recovery_rate, 1.0/time_steps_in_day);
		enum class State {
			Infected, Recovered, Susceptible
		};

		vector<vector<int>> infected_to_susceptible(n);
		vector<State> state(n);
		vector<int> infected;

		int nsusceptible=0, ninfected=0, nrecovered=0;
		for (int i=0; i<n; i++) {
			if (uniform_real_distribution<>()(rng) < initial_infection_prob) {
				state[i]=State::Infected;
				infected.push_back(i);
				ninfected++;
			} else {
				state[i]=State::Susceptible;
				nsusceptible++;
			}
		}

		vector<int> degs(n,0);
		auto process_infected_to_susceptible = [&](int i) {
			if (stay_home) {
				int k = uniform_int_distribution<>(1,K)(rng);
				for (int x: edges[i]) degs[x]=1;

				set<array<int,2>> left;
				for (int x: edges[i]) {
					for (int y: edges[x]) {
						if (degs[y]>0) degs[x]++;
					}

					left.insert({ degs[x],x });
				}

				while (left.size() && (*left.begin())[0]<k) {
					int x = (*left.begin())[1];
					left.erase(left.begin());
					degs[x]=0;

					for (int y: edges[x]) {
						if (degs[y]>0) {
							left.erase({ degs[y], y });
							left.insert({ --degs[y], y });
						}
					}
				}

				for (int x: edges[i]) degs[x]=0;

				for (auto [deg, x]: left) {
					if (state[x]==State::Susceptible)
						infected_to_susceptible[i].push_back(x);
				}
			} else {
				for (int x: edges[i]) {
					if (state[x]==State::Susceptible)
						infected_to_susceptible[i].push_back(x);
				}
			}

			shuffle(infected_to_susceptible[i].begin(), infected_to_susceptible[i].end(), rng);
		};

		for (int i=0; i<n; i++) {
			if (state[i]==State::Infected) process_infected_to_susceptible(i);
		}

		vector<int> new_infected;
		int plot_i=0;
		for (int day=1; day<=num_days; day++) {
			for (int j=0; j<time_steps_in_day; j++) {
				if (j%plot_every_time_step==0) {
					total_infected[plot_i] += ninfected;
					total_susceptible[plot_i] += nsusceptible;
					plot_i++;
				}

				new_infected.clear();

				for (int k: infected) {
					int num_to_infect = binomial_distribution<>(infected_to_susceptible[k].size(), time_step_infection_rate)(rng);

					while (num_to_infect--) {
						int y = infected_to_susceptible[k].back();
						infected_to_susceptible[k].pop_back();

						state[y]=State::Infected;
						new_infected.push_back(y);
						nsusceptible--, ninfected++;

						for (int z: edges[y]) {
							if (z!=k && state[z]==State::Infected) {
								auto it = find(infected_to_susceptible[z].begin(), infected_to_susceptible[z].end(), y);
								if (it!=infected_to_susceptible[z].end())
									infected_to_susceptible[z].erase(it);
							}
						}

						process_infected_to_susceptible(y);
					}

					if (uniform_real_distribution<>()(rng)<time_step_recovery_rate) {
						state[k]=State::Recovered;
						nrecovered++, ninfected--;
					} else {
						new_infected.push_back(k);
					}
				}

				new_infected.swap(infected);
			}

			cout<<"Trial "<<trial<<", day "<<day<<": "<<nsusceptible<<" susceptible, "<<ninfected<<" infected, "<<nrecovered<<" recovered\n";
		}
	}

	ofstream plot("./plot.gp");
	plot<<"$DATA<<e\n";

	int plot_i=0;
	double fac = 100.0/n/num_trials;
	for (int day=1; day<=num_days; day++) {
		for (int j=time_steps_in_day*day; j<time_steps_in_day*(day+1); j+=plot_every_time_step) {
			plot<<1.0*j/time_steps_in_day
				<<" "<<total_infected[plot_i]*fac
				<<" "<<total_susceptible[plot_i]*fac
				<<" "<<100 - (total_infected[plot_i] + total_susceptible[plot_i])*fac
				<<"\n";
			plot_i++;
		}
	}

	plot<<"e\n";

	plot<<format(R"EOF(
set terminal pdfcairo lw 3.0
set border lw 0.2
set output "plot.pdf"
set term pdfcairo font "Playfair Display bold,12"
set key inside center right

set title "Graph of small-world-based SIR model, averaged over 10 trials{2}"
set xlabel "Day"
set ylabel "% of population (n={1})"
set yrange [0:100]
set xrange [1:{0}]

set xtics 1,5 font ",8"
set ytics format "%.0f%%"

set grid ls 2 linecolor "gray90" linewidth 0.5 dashtype solid

plot $DATA using 1:2 with lines title 'SIR model (infected)', $DATA using 1:3 with lines title 'SIR model (susceptible)', $DATA using 1:4 with lines title 'SIR model (recovered)'
	)EOF", num_days, n, stay_home ? "\\n{/*0.8 (with some students staying home)}" : "");
}
 