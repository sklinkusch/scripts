# Mehrheit

## Perl

The Perl script is used with
```shell
./mehrheit.pl A B PartyA 221 PartyB 125 PartyC 76 PartyD 32
```
### Parameters
#### Mode A
- **A** more than 50%, no excluded coalitions
- **B** more than 50%, with excluded coalitions
- **C** more than 40%, no excluded coalitions
- **D** more than 40%, with excluded coalitions
- **E** no excluded coalitions
- **F** with excluded coalitions

#### Mode B
- **A** no oversize coalitions
- **B** oversize coalitions

#### Parties and seats
After the two modes the party names and the respective number of seats follow, separated by spaces.

#### Excluded combinations
Excluded combinations of two parties are written in one line (separated by a space) in the file `exclusion.list`. 

## Typescript

```shell
npx ts-node ./mehrheit.ts [options] [optionValues]
```

### Options
#### Options without Values
- `--size40`: more than 40%
- `--sizeAll`: all combinations
- `--oversize`: also oversize coalitions

#### Options with values
- `--exclude EXCLUDE-LIST`: exclude combinations given in `EXCLUDE-LIST`. `EXCLUDE-LIST`  is a JSON array of arrays, encapsulated in single quotes. Each subarray contains the parties that will not form a coalition together
- `--values VALUES`: gives the parties with the respective seats, given in the JSON object `VALUES`, encapsulated in single quotes.