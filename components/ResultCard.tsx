import React from 'react';
import { LucideIcon } from 'lucide-react';

interface ResultCardProps {
  label: string;
  value: string | number;
  unit?: string;
  icon: LucideIcon;
  subValue?: string;
}

export const ResultCard: React.FC<ResultCardProps> = ({ label, value, unit, icon: Icon, subValue }) => {
  return (
    <div className="bg-slate-100 p-4 rounded-xl border border-slate-200 flex items-center space-x-4 shadow-sm">
      <div className="p-3 bg-white rounded-lg shadow-sm border border-slate-200">
        <Icon className="w-6 h-6 text-indigo-600" />
      </div>
      <div>
        <p className="text-xs font-bold text-slate-500 uppercase">{label}</p>
        <p className="text-2xl font-bold text-slate-900 tracking-tight">
          {value} <span className="text-sm text-slate-500 font-normal">{unit}</span>
        </p>
        {subValue && <p className="text-xs text-slate-500 mt-1">{subValue}</p>}
      </div>
    </div>
  );
};