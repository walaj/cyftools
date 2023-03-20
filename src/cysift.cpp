#include "cell_column.h"
#include "cell_table.h"

int main() {
  
  CellTable table("-", "markers_orion.json", true);

  table.PrintPearson(false);
  
  /*  std::thread producer_thread(&CellTable::producer, &table, "-");
  std::thread consumer_thread(&CellTable::consumer, &table);
  
  producer_thread.join();
  consumer_thread.join();
  */
  //std::cerr << table;
   
  return 0;
}
