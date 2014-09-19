import std.stdio, std.algorithm, std.range;

void main(string[] args) {

    bool keep = true;

    auto input = File(args[1]);
    auto output = File(args[1] ~ ".filtered", "w+");
    foreach (line; input.byLine()) {
      if (line.length > 0 && line[0] == '#') {
        auto fn = line.splitter(' ').drop(2).front;
        keep = !fn.startsWith(`"/usr/`);
      } else {
        if (keep)
          output.writeln(line);
      }
    }
}
