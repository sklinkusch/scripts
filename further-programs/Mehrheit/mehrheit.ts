const args = process.argv.slice(2);

// function that returns all possible combinations of an existing array
function getCombinations(array: string[]): (string[])[] {
  let combi = [];
  let temp = [];
  const slent = Math.pow(2, array.length);
  for (let i = 0; i < slent; i++) {
    temp = [];
    for (let j = 0; j < array.length; j++) {
      if (i & Math.pow(2, j)) {
        temp.push(array[j]);
      }
    }
    if (temp.length > 0) {
      combi.push(temp);
    }
  }
  combi.sort((a, b) => a.length - b.length);
  return combi;
}

let allowOverSize: boolean = false;
let sizeAll: boolean = false;
let size40: boolean = false;
let size50: boolean = true;
let useExclusionList: boolean = false;
let exclusionList: string[][] = [];
let values: { [key: string]: number } = {};

if (args.includes('--oversize')) allowOverSize = true;
if (args.includes('--all')) sizeAll = true;
if (args.includes('--40')) size40 = true;
if (args.includes('--exclude')) {
  useExclusionList = true;
  const excludeKeywordIndex = args.indexOf('--exclude');
  const exclusionListIndex = excludeKeywordIndex + 1;
  const exclusionListJson = args[exclusionListIndex];
  exclusionList = JSON.parse(exclusionListJson);
}
if (args.includes('--values')) {
  const valueKeywordIndex = args.indexOf('--values');
  const valueListIndex = valueKeywordIndex + 1;
  const valueListJson = args[valueListIndex];
  values = JSON.parse(valueListJson);
}
// sort parties by size
const parties = Object.keys(values);
const seats = Object.values(values);
const sortedParties = parties.sort((a, b) => {
  return values[b] - values[a];
});
// calculate all seats and majority
const allSeats = seats.reduce((sum, acc) => sum + acc, 0);
const majority = allSeats % 2 === 0 
  ? allSeats / 2 + 1
  : allSeats / 2 + 0.5;
const minSize40 = 0.4 * allSeats;

// find all combinations
const combinations = getCombinations(sortedParties);
// exclude combinations that do not have sufficient votes
const minimumVotes = sizeAll
  ? 0
  : size40
    ? minSize40
    : majority;
let validCombinations = combinations.filter(combination => {
  let totalVotes = combination.reduce((sum, acc) => sum + values[acc], 0);
  return totalVotes >= minimumVotes;
});
// exclude combinations that are on the exclusion list
if (useExclusionList) {
  let k = 0;
  while (k < validCombinations.length) {
    const currentCombination = validCombinations[k];
    const isExcluded = exclusionList.reduce((acc, comb) => {
      if (acc) return true;
      if (comb.every(party => currentCombination.includes(party))) return true;
      return false;
    }, false);
    if (isExcluded) {
      validCombinations.splice(k, 1);
    } else {
      k++;
    }
  }
}
// exclude oversize coalitions if allowOverSize is not true
if (allowOverSize) {
  // do nothing
} else {
  validCombinations.forEach((comb, index) => {
    let i = index + 1;
    while(i < validCombinations.length) {
      const allPartiesIncluded = comb.reduce((acc, party) => {
        return acc && validCombinations[i].includes(party);
      }, true);
      if (allPartiesIncluded) {
        validCombinations.splice(i, 1);
      } else {
        i++;
      }
    }
  });
}

const validCombinationArray = validCombinations.map(combination => {
  return {
    name: combination.join('+'),
    total: combination.reduce((sum, party) => sum + values[party], 0),
    percent: (100 * combination.reduce((sum, party) => sum + values[party], 0) / allSeats).toFixed(1) + '%',
  };
});
const sortedValidCombinations = validCombinationArray.sort((a, b) => a.total - b.total);
console.table(sortedValidCombinations);